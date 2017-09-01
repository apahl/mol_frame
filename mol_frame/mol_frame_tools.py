#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############
MolFrame Tools
##############

*Created on Thu Aug 31, 2017 by A. Pahl*

Tools for working with MolFrame dataframes.
"""

import time
import os.path as op
# from collections import Counter

import pandas as pd
# import numpy as np

from rdkit.Chem import AllChem as Chem
# from rdkit.Chem import Draw

from IPython.core.display import HTML

from .viewers import df_html, view

try:
    from misc_tools import apl_tools
    AP_TOOLS = True
    print("* reimported mol_frame_tools")
    #: Library version
    # VERSION = apl_tools.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    # print("{:45s} (commit: {})".format(__name__, VERSION))

except ImportError:
    AP_TOOLS = False
    print("{:45s} ({})".format(__name__, time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__)))))

try:
    # Try to import Avalon so it can be used for generation of 2d coordinates.
    from rdkit.Avalon import pyAvalonTools as pyAv
    USE_AVALON_2D = True
except ImportError:
    print("  * Avalon not available. Using RDKit for 2d coordinate generation.")
    USE_AVALON_2D = False

DEBUG = True


def debug_print(txt, val):
    if DEBUG:
        txt = txt + ":"
        print("DEBUG   {:20s}".format(txt), val)


class MolFrame():
    def __init__(self, init=None, log=True):
        if init is None:
            self.data = pd.DataFrame()
        else:
            self.data = pd.DataFrame(init)
        self.has_mols = False


    def __getitem__(self, item):
        res = self.data[item]
        if isinstance(res, pd.DataFrame):
            result = MolFrame()
            result.data = res
            result.print_log("subset")
        else:
            result = res
        return result


    def __getattr__(self, name):
        """Try to call undefined methods on the underlying pandas DataFrame."""
        def method(*args, **kwargs):
            res = getattr(self.data, name)(*args, **kwargs)
            if isinstance(res, pd.DataFrame):
                result = MolFrame()
                result.data = res
                result.print_log(name)
            else:
                result = res
            return result
        return method


    def print_log(self, component, add_info=""):
        if self.log:
            print_log(self.data, component, add_info)


    def show(self, include_smiles=False, drop=[], **kwargs):
        smiles_col = kwargs.get("smiles_col", "Smiles")
        if not include_smiles:
            drop.append(smiles_col)
        if len(drop) > 0:
            self.drop_cols(drop)
        return HTML(df_html(self.data))


    def view(self, drop=[], keep=[], fn="molframe.html", **kwargs):
        self.add_mols()
        return view(self.data, drop=drop, keep=keep, fn=fn, **kwargs)


    def head(self, n=5):
        res = self.data.head(n)
        result = MolFrame()
        result.data = res
        result.print_log("head")
        return result


    def drop_cols(self, cols, inplace=False):
        """Drops the list of columns from the DataFrame.
        Listed columns that are not present in the DataFrame are simply ignored
        (no error is thrown)."""
        if inplace:
            drop_cols(self.data, cols, inplace=True)
            self.print_log("drop cols (inplace)")
        else:
            result = MolFrame()
            result.data = drop_cols(self.data, cols, inplace=False)
            result.print_log("drop cols")
            return result


    def keep_cols(self, cols, inplace=False):
        if inplace:
            self.data = self.data[cols]
            self.print_log("keep cols (inplace)")
        else:
            result = MolFrame()
            result.data = self.data[cols]
            result.print_log("keep cols")
            return result


    def load(self, fn, sep="\t"):
        """Read one or multiple result files and concatenate them into one MolFrame.
        `fn` is a single filename (string) or a list of filenames."""
        self.data = load(fn, sep=sep).data
        self.print_log("load data")


    def write_csv(self, fn, parameters=None, sep="\t"):
        result = self.data.copy()
        if isinstance(parameters, list):
            result = result[parameters]
        result.to_csv(fn, sep=sep, index=False)


    def write_pkl(self, fn):
        self.data.to_pickle(fn)


    def add_mols(self, smiles_col="Smiles", force=False):
        def _mol_from_smiles(smi):
            return pd.Series(mol_from_smiles(smi))
        if force or not self.has_mols:
            self.data["Mol"] = self.data[smiles_col].apply(_mol_from_smiles)
            self.has_mols = True


def check_2d_coords(mol, force=False):
    """Check if a mol has 2D coordinates and if not, calculate them."""
    if not force:
        try:
            mol.GetConformer()
        except ValueError:
            force = True  # no 2D coords... calculate them

    if force:
        if USE_AVALON_2D:
            pyAv.Generate2DCoords(mol)
        else:
            mol.Compute2DCoords()


def print_log(df, component, add_info=""):
    component = component + ":"
    if len(add_info) > 0:
        add_info = "    ({})".format(add_info)
    print("* {:22s} ({:5d} | {:4d}){}".format(component, df.shape[0], df.shape[1], add_info))


def load(fn, sep="\t"):
    """Read one or multiple result files and concatenate them into one MolFrame.
    `fn` is a single filename (string) or a list of filenames."""
    result = MolFrame()
    if isinstance(fn, list):
        result.data = pd.concat((pd.read_csv(f, sep=sep) for f in fn))
    else:
        result.data = pd.read_csv(fn, sep=sep)

    result.print_log("load molframe")
    return result


def mol_from_smiles(smi, calc_2d=True):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles("*")
    else:
        if calc_2d:
            check_2d_coords(mol)
    return mol


def drop_cols(df, cols, inplace=False):
    """Drops the list of columns from the DataFrame.
    Listed columns that are not present in the DataFrame are simply ignored
    (no error is thrown)."""
    df_keys = df.keys()
    drop = [k for k in cols if k in df_keys]
    if inplace:
        df.drop(drop, axis=1, inplace=True)
    else:
        result = df.drop(drop, axis=1)
        return result
