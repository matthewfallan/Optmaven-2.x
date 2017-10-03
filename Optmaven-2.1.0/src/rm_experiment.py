import sys

import standards

experiment = sys.argv[1]
standards.safe_rm_experiment(experiment)
