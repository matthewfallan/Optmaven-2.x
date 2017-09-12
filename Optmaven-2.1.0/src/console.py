""" This module controls output to the console. """

import os

import standards


def disp(text):
    print(wrap(text))


def clear():
    os.system("clear")


def wrap(text, width=None):
    # FIXME: wrap text to standards.ConsoleWidth lines.
    if width is None:
        width = standards.ConsoleWidth
    return text


def disp_list(list_, sep=", ", max_display=standards.MaxListDisplay, get_string=False):
    if max_display is None:
        max_display = len(list_)
    max_display = min(max_display, len(list_))
    if max_display < len(list_):
        items = list_[: max_display - 1] + ["...", list_[-1]]
    else:
        items = list_
    joined = sep.join(map(str, items))
    if get_string:
        return wrap(joined)
    else:
        disp(joined)
