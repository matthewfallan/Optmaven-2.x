
import subprocess

def submit(file_name, command, walltime):
    secs = int(walltime % 60)
    walltime_mins = int((walltime - secs) / 60)
    mins = walltime_mins % 60
    walltime_hours = (walltime_mins - mins) / 60
    hours = walltime_hours
    walltime_text = "{0:>2}:{0:>2}:{0:>2}".format(hours, mins, secs)
    script = """#!/bin/sh
#PBS -q lionxf
#PBS -l walltime={t}
{c}""".format(t=walltime_text, c=command)
    with open(file_name, "w") as f:
        f.write(script)
    subprocess.call(["chmod", "u+x", file_name])
    status = subprocess.call([file_name])
    #status = subprocess.call(["qsub", file_name])
