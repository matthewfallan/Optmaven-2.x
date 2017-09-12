__doc__ = """ This module provides a consistent interface to the MAPs
database. """

import cPickle as pkl
import itertools
import os

import numpy as np

import standards

coords_files = dict()
gaps = [8]
for gap in gaps:
    coords_files[gap] = os.path.join(standards.MapsDirectory, "all_coords_{}.csv".format(
            gap))

#FIXME
chains = standards.MapsChains
cdrs = standards.MapsCdrs
categories = ["{}{}".format(chain, cdr) for chain, cdr in itertools.product(chains, cdrs)]
category_directories = [os.path.join(standards.MapsDirectory, category) for category in categories]
parts = {os.path.splitext(part)[0]: os.path.join(category, part) for category in category_directories for part in os.listdir(category)}
#part_category = {part: category for category, category_parts in parts.iteritems()
#        for part in category_parts}
#coords = dict()
#for gap in gaps:
#    coords[gap] = {line.split(",")[0]: np.array(map(float, line.split(",")[1:]))
#            for line in open(coords_files[gap])}


def get_part_category(part):
    return part_category[part]


def get_coordinates(part, gap):
    # FIXME: add support for different gap penalties in the future
    c = coords[gap].get(part)
    # FIXME: this is just a temporary hack for duplicate MAPs parts
    if c is None:
        c = coords[gap][part.split("_")[0] + "_" + str(int(part.split("_")[1]) + 1)]
    return c


def get_parts():
    """ Return {category1: [part1, part2, ... ], category2: [ ... ] ... } """
    return {category: sorted(parts.keys(), key=lambda x: int(x.split("_")[1]))
            for category, parts in parts.iteritems()}


def list_categories():
    return sorted(parts.keys())


def list_parts(categories=None):
    """ Return a flat list of the names of MAPs parts. A category or
    multiple categories may be speccified. """
    if categories is None:
        categories = list_categories()
    elif isinstance(categories, str):
        categories = [categories]
    parts_list = [part for category in categories for part in sorted(parts[
            category])]
    return parts_list


def get_part_file(part):
    """ Return the path to the file for the part. """
    return parts[get_part_category(part)][part]
