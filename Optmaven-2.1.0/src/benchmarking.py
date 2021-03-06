from collections import OrderedDict
import datetime
import os
import subprocess

import standards


class BlankTimeFileError(IOError):
    pass


class Task(object):
    def __init__(self, purpose, detail=None, time_stamp=None):
        self.purpose = purpose
        self.detail = detail
        if time_stamp is None:
            self.time_stamp = datetime.datetime.now()
        else:
            self.time_stamp = time_stamp

    def get_time_stamp(self):
        return self.time_stamp


class Time(Task):
    def __init__(self, purpose, time, detail=None, time_stamp=None):
        Task.__init__(self, purpose, detail, time_stamp)
        self.time = time

    def to_dict(self):
        info = {
            "Type": "Time",
            "Purpose": self.purpose,
            "Detail": self.detail,
            "Time Stamp": self.time_stamp.strftime(standards.DatetimeFormat)
        }
        for _type in standards.UnixTimeCodes:
            info[_type] = self.time[_type]
        return info
            

class TimeFile(Task):
    def __init__(self, _file, purpose, detail=None):
        Task.__init__(self, purpose, detail)
        self.file = _file

    def to_dict(self):
        return Time(self.purpose, parse_time_file(self.file), self.detail, self.time_stamp).to_dict()
            

class DriveUsage(Task):
    def __init__(self, experiment, detail=None):
        Task.__init__(self, experiment.purpose, detail)
        du = subprocess.Popen(["du", "-s", experiment.directory], stdout=subprocess.PIPE)
        out, error = du.communicate()
        drive_usage, directory = out.split("\t")
        self.drive_usage = float(drive_usage)
        
    def to_dict(self):
        info = {
            "Type": "Drive Usage",
            "Purpose": self.purpose,
            "Detail": self.detail,
            "Drive Usage": self.drive_usage,
            "Time Stamp": self.time_stamp.strftime(standards.DatetimeFormat)
        }
        return info
 

def parse_time_file(_file):
    with open(_file) as f:
        lines = f.readlines()
    if len(lines) == 0:
        raise BlankTimeFileError("Found blank time file: {}".format(_file))
    if len(lines) != len(standards.UnixTimeCodes):
        raise ValueError("Need {} lines in time file {}".format(len(standards.UnixTimeCodes), _file))
    return {time_type: float(line) for line, time_type in zip(lines, standards.UnixTimeCodes)}
