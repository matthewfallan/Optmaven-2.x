""" This module provides functions for user input. """

import copy
import os
import string

from Bio.PDB import PDBList

from console import clear, disp, disp_list, wrap
import standards


def get(prompt):
    return raw_input(wrap(prompt))

def get_number(prompt, min_value=None, max_value=None, value_type=float):
    if value_type not in (float, int):
        raise ValueError("Numeric type must be float or int, not {}.".format(value_type))
    do = True
    while do:
        try:
            answer = value_type(get(prompt))
        except ValueError:
            pass
        else:
            if min_value is not None and answer < min_value:
                disp("Value must be greater than or equal to {}.".format(min_value))
            elif max_value is not None and answer > max_value:
                disp("Value must be less than or equal to {}.".format(max_value))
            else:
                do = False
    return answer
            

def get_file(prompt, directory, new_file=False, fetch_pdb=False):
    if not os.path.isdir(directory):
        raise OSError("Directory does not exist: {}".format(directory))
    if new_file:
        do = True
        while do:
            base_file = get(prompt)
            full_file = os.path.join(directory, base_file)
            if os.path.exists(full_file):
                disp("Cannot create file, path already exists: {}".format(full_file))
            else:
                if standards.is_path_component(base_file):
                    file_directory = os.path.dirname(full_file)
                    if os.path.isdir(file_directory):
                        do = False
                    else:
                        disp("Cannot create file, directory does not exist: {}".format(file_directory))
                else:
                    disp("Only the following characters are permitted in path components: {}".format(standards.AllowedPathCharacters))
        return full_file
    else:
        if fetch_pdb:
            fetch_pdir = directory
        else:
            fetch_pdir = None
        return os.path.join(directory, select_one_from_list(prompt, os.listdir(directory), fetch_pdir=fetch_pdir))


def get_files(prompt, directory, min_number, max_number):
    return [os.path.join(directory, file_) for file_ in select_from_list(prompt, os.listdir(directory), min_number, max_number)]


def get_yn(prompt):
    yes = ["yes", "y"]
    no = ["no", "n"]
    answer = None
    while answer not in yes + no:
        answer = get(prompt).lower()
    return answer in yes


def parse_numeric_selection(selection_text):
    """ Return a set of whole numbers based on a selection string. The string should be formatted with commas and dashes, e.g. "1,3,7-9" returns [1, 3, 7, 8, 9]. Spaces do not matter. """
    selection_items = "".join(selection_text.split()).split(",")
    selection_numbers = list()
    for item in selection_items:
        if "-" in item:
            bounds = item.split("-")
            if len(bounds) != 2:
                raise ValueError('Cannot parse range: {}. Must be formatted "x-y".'.format(item))
            selection_numbers.extend(range(int(bounds[0]), int(bounds[1]) + 1))
        else:
            selection_numbers.append(int(item))
    return selection_numbers


def select_from_list(prompt, list_, min_number=None, max_number=None, names=None, fetch_pdir=None):
    n = len(list_)
    if isinstance(list_, dict):
        if names is None:
	    names = map(str, list_.keys())
            list_ = list_.values()
        else:
            list_ = list_.keys()
    if names is None:
        names = map(str, list_)
    else:
        if len(names) != n:
            raise ValueError("Length of names must equal length of list.")
    names = map(string.strip, names)
    if min_number is None:
        min_number = 0
    else:
        min_number = max(int(min_number), 0)
    if max_number is None:
        max_number = n
    else:
        max_number = min(int(max_number), n)
    if min_number > max_number:
        raise IOError("Minimum selction {} cannot be greater than maximum selection {}.".format(min_number, max_number))
    if min_number > n:
        raise IOError("Cannot select at minimum {} items from list with length {}.".format(min_number, n))
    if min_number == n:
        return list_
    items = None
    permitted = "Options include the following:\n{}".format(disp_list(names, get_string=True))
    permitted_displayed = False
    while items is None:
        selection_string = get(prompt)
        if selection_string.strip() == "":
            disp(permitted)
            permitted_displayed = True
        else:
            selection_keys = split_selection_string(selection_string)
            try:
                indexes = [index for key in selection_keys for index in get_selected_indexes(names, key)]
            except ValueError:
                # If selecting a PDB, try to import it.
                if fetch_pdir is not None:
                    disp("Trying to fetch missing PDB files.")
                    successful_fetch = True
                    for key in selection_keys:
                        pdb_code = os.path.splitext(key)[0].upper()
                        #disp(pdb_code)
                        PDBList().retrieve_pdb_file(pdb_code, pdir=fetch_pdir)
                        if os.path.isfile(os.path.join(fetch_pdir, key)):
                            list_.append(key)
                            names.append(str(key))
                            n += 1
                        else:
                            successful_fetch = False
                    if successful_fetch:
                        try:
                            indexes = [index for key in selection_keys for index in get_selected_indexes(names, key)]
                        except ValueError:
                            select_items = False
                        else:
                            select_items = True
                    else:
                        select_items = False
                else:
                    select_items = False
            else:
                select_items = True
            if select_items:
                selected_items = [list_[index] for index in indexes]
                if len(selected_items) >= min_number and len(selected_items) <= max_number:
                    items = selected_items
                else:
                    disp("You must select from {} to {} items.".format(min_number, max_number))
            else:
                disp('There was a problem selecting the following items: "{}"'.format(selection_string))
                if not permitted_displayed:
                    disp(permitted)
                    permitted_displayed = True
 


    return items


def split_selection_string(string, sep=","):
    return [item.strip() for item in string.split(sep)]


def get_selected_indexes(list_, key, range_sep="-"):
    if key == standards.SelectAll:
        return range(len(list_))
    elif key == standards.SelectNone:
        return list()
    elif key.count(range_sep) == 1:
        # If the selection key indicates a range, then use the given items as bounds.
        lower, upper = [unescape_item(item.strip()) for item in key.split(range_sep)]
        # Return all items between them.
        i_lower = list_.index(lower)
        i_upper = list_.index(upper)
        return range(i_lower, i_upper + 1) 
    else:
        # Otherwise, just use the one index.
        return [list_.index(unescape_item(key))]


def unescape_item(item):
    if item.startswith(standards.EscapeCharacter):
        return item[1:]
    else:
        return item


def select_one_from_list(prompt, list_, names=None, fetch_pdir=None):
    return select_from_list(prompt, list_, 1, 1, names, fetch_pdir)[0]


if __name__ == "__main__":
    get_files("files", os.getcwd(), 2, 4)
