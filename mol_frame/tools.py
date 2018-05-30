#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#####
Tools
#####

*Created on Mon, Jan 22 2018, 12:00 by A. Pahl*

Helper Tools, e.g. for loading the configuration.
"""

import os
import os.path as op
import yaml
import math

import numpy as np


def load_config(conf="config"):
    """Loads configuration from default location and
    returns config object.
    Known configuration files are `config.yaml` and `plates.yaml`.
    Raises error when the file could not be loaded."""
    assert conf in ["config", "plates"]

    if "HOME" in os.environ:
        conf_fn = op.join(os.environ["HOME"], ".config",
                          "mol_frame", "{}.yaml".format(conf))
    elif "HOMEPATH" in os.environ:  # Windows
        conf_fn = op.join(os.environ["HOMEPATH"],
                          "mol_frame", "{}.yaml".format(conf))
    try:
        with open(conf_fn, 'r') as ymlfile:
            config = yaml.load(ymlfile)
    except FileNotFoundError:
        print("Configuration file {}.yaml not found.".format(conf))
        print("`load_resources()` will not work.")
        print("Have a look at the *.yaml files in the `conf` folder")
        print("for templates and locations.")
        config = {}
    return config


def unit_factor(unit):
    """Return the factor corresponding to the unit, e.g. 1E-9 for nM.
    Known units are: mM, uM, nM, pM. Raises ValueError for unknown unit."""
    units = ["mm", "um", "nm", "pm"]
    pos = units.index(unit.lower()) + 1
    factor = 10 ** -(pos * 3)
    return factor


def pic50(ic50, unit=None, digits=3):
    """Calculate pIC50 from IC50. Optionally, a unit for the input IC50 value may be given.
    Known units are: mM, uM, nM, pM"""
    if unit is not None:
        ic50 *= unit_factor(unit)
    return np.round(-math.log10(ic50), decimals=digits)


def ic50(pic50, unit=None, digits=3):
    """Calculate IC50 from pIC50. Optionally, a unit for the returned IC50 value may be given.
    Known units are: mM, uM, nM, pM"""
    ic50 = 10 ** (-pic50)
    if unit is not None:
        ic50 /= unit_factor(unit)
    return np.round(ic50, digits)
