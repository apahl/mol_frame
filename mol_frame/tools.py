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
