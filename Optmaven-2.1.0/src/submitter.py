import cPickle as pkl
import os
import re
import subprocess
import sys
import tempfile

import standards


class PbsBatchSubmitter(object):
    def __init__(self, experiment):
        self.experiment = experiment
        self.directory = tempfile.mkdtemp(dir=self.experiment.get_temp(), prefix="pbs_")
        self.file = os.path.join(self.directory, "pbs.pickle")
        self.jobs_file_prefix = os.path.join(self.directory, standards.PbsJobFilePrefix)
        self.jobs = dict()

    def save(self):
        with open(self.file, "w") as f:
            pkl.dump(self, f)

    def get_jobs_file(self, array):
        return os.path.join(self.directory, "{}{}".format(standards.PbsJobFilePrefix, array))

    def get_jobs_files(self):
        return map(self.get_jobs_file, range(1, self.array + 1))

    def collect_garbage(self):
        self.experiment.safe_rmtree(self.directory)

    def submit(self, program, args, jobs):
        # jobs is a dictionary where the keys are the arguments that must be passed to experiment.py and the values are the files that are generated when the jobs have completed.
        self.program = program
        self.args = args
        self.jobs.update(jobs)
        # Remove all previous jobs files before submitting new jobs.
        self.collect_garbage()
        try:
            os.mkdir(self.directory)
        except OSError:
            pass
        unfinished = [job for job, _file in self.jobs.iteritems() if not os.path.isfile(_file)]
        if len(unfinished) > 0:
            self.array = 0
            collection = list()
            for i, job in enumerate(unfinished):
                collection.append(job)
                if len(collection) >= self.experiment.batch_size or i == len(unfinished) - 1:
                    self.array += 1
                    jobs_file = self.get_jobs_file(self.array)
                    with open(jobs_file, "w") as f:
                        f.write(" ".join(map(str, collection)))
                    collection = list()
            command = "{} {} {}{}".format(program, " ".join(map(str, args)), self.jobs_file_prefix, standards.PbsArrayId)
            if command.startswith("rm"):
                raise OSError("Submitting 'rm' is forbidden.")
            handle, self.script_file = tempfile.mkstemp(dir=self.directory, prefix="script_", suffix=".sh")
            os.close(handle)
            write_script(self.script_file, command, self.experiment.walltime, self.array)
            call = [standards.PbsQsub, self.script_file]
            p = subprocess.Popen(call, stdout=subprocess.PIPE)
            stdout, stderr = p.communicate()
            job_id_match = re.match("[0-9]+\[\]", stdout)
            if job_id_match is None:
                raise OSError("Unable to implement callback for job id {}".format(stdout))
            job_id = job_id_match.group()
            command = "{} {} {}".format(standards.PythonCommand, os.path.realpath(__file__), self.file)
            handle, self.callback_file = tempfile.mkstemp(dir=self.directory, prefix="callback_", suffix=".sh")
            os.close(handle)
            self.save()
            submit(self.callback_file, command, self.experiment.walltime, options={"-W": "depend=afteranyarray:{}".format(job_id)})
        else:
            self.collect_garbage()
            self.experiment.run_next()


def script_initial(queue, walltime):
    secs = int(walltime % 60)
    walltime_mins = int((walltime - secs) / 60)
    mins = walltime_mins % 60
    walltime_hours = (walltime_mins - mins) / 60
    hours = walltime_hours
    walltime_text = "{0:>2}:{0:>2}:{0:>2}".format(hours, mins, secs)
    lines = ["#!/bin/sh",
             "#PBS -q {}".format(queue),
             "#PBS -l walltime={}".format(walltime_text)]
    return lines


def write_script(file_name, command, walltime, array=0):
    lines = script_initial(standards.PbsQueue, walltime)
    if array > 0:
        lines.append("#PBS -t {}-{}".format(1, array))
    if isinstance(command, list):
        lines.extend(command)
    else:
        lines.append(command)
    script = "\n".join(map(str, lines))
    with open(file_name, "w") as f:
        f.write(script)
    subprocess.call(["chmod", "u+x", file_name])


def submit(file_name, command, walltime, options=None, queue=True):
    write_script(file_name, command, walltime)
    if options is None:
        options = dict()
    if queue:
        call = [standards.PbsQsub]
    else:
        call = list()
    [call.extend([k, v]) for k, v in options.items()]
    call.append(file_name)
    status = subprocess.call(call)


if __name__ == "__main__":
    sub_file = sys.argv[1]
    with open(sub_file) as f:
        sub = pkl.load(f)
    sub.submit(sub.program, sub.args, dict())
