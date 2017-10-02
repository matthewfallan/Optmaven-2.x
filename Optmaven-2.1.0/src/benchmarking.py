from collections import OrderedDict

import standards


class Task(object):
    def __init__(self, purpose, detail=None):
        self.purpose = purpose
        self.detail = detail


class Time(Task):
    def __init__(self, purpose, time, detail=None):
        Task.__init__(self, purpose, detail)
        self.time = time

    def to_dict(self):
        info = OrderedDict([
            ("Type", "Time"),
            ("Purpose", self.purpose),
            ("Detail", self.detail)
        ])
        for _type in standards.UnixTimeCodes:
            info[_type] = self.time[_type]
        return info
            

class TimeFile(Task):
    def __init__(self, _file, purpose, detail=None):
        Task.__init__(self, purpose, detail)
        self.file = _file

    def to_dict(self):
        with open(self.file) as f:
            times = {time_type: float(line) for line, time_type in zip(f, standards.UnixTimeCodes)}
        return Time(self.purpose, times, self.detail).to_dict()    
            
