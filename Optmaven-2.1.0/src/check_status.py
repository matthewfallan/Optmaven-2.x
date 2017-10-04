import cPickle as pkl
import os

import console
import standards


console.disp("{:<39} {}".format("Experiment", "Status"))
for experiment in os.listdir(standards.ExperimentsDirectory):
    directory = os.path.join(standards.ExperimentsDirectory, experiment)
    errors = os.path.join(directory, "errors.txt")
    pickle = os.path.join(directory, ".temp", "{}.pickle".format(experiment))
    summary = os.path.join(directory, "Summary.txt")
    results = os.path.join(directory, "Results.csv")
    if os.path.isfile(errors):
        status = "ERROR"
    else:
        try:
            with open(pickle) as f:
                exp = pkl.load(f)
        except IOError:
            if all([os.path.isfile(x) for x in [summary, results]]):
                status = "Completed"
            else:
                raise IOError("Cannot find results and summary or pickle file for Experiment {}".format(experiment))
        else:
            task, status = exp.get_task(exp.status)
    console.disp("{:<39} {}".format(experiment, status))