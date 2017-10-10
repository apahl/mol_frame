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
import pickle  # for b64 encoded mol objects
import base64 as b64
import os.path as op
# from collections import Counter

import pandas as pd
# import numpy as np

from rdkit.Chem import AllChem as Chem
# from rdkit.Chem import Draw

from .viewers import df_html, view

try:
    from . import resource_paths as cprp, nb_tools as nbt
except ImportError:
    from . import resource_paths_templ as cprp
    print("* Resource paths not found, stub loaded.")
    print("  Automatic loading of resources will not work,")
    print("  please have a look at resource_paths_templ.py")

try:
    from misc_tools import apl_tools
    AP_TOOLS = True
    print("* reimported mol_frame")
    #: Library version
    VERSION = apl_tools.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} (commit: {})".format(__name__, VERSION))

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


IPYTHON = nbt.is_interactive_ipython()

if IPYTHON:
    from IPython.core.display import HTML

DEBUG = False


def debug_print(txt, val):
    if DEBUG:
        txt = txt + ":"
        print("DEBUG   {:20s}".format(txt), val)


class MolFrame(object):
    def __init__(self, init=None, log=True):
        if init is None:
            self.data = pd.DataFrame()
        else:
            self.data = pd.DataFrame(init)
        self.inplace = False
        self.has_mols = False
        self.id_col = "Compound_Id"
        self.smiles_col = "Smiles"
        self.mol_col = "Mol"
        self.b64_col = "Mol_b64"


    def _pass_properties(self, mol_frame):
        mol_frame.inplace = self.inplace
        mol_frame.has_mols = self.has_mols
        mol_frame.id_col = self.id_col
        mol_frame.smiles_col = self.smiles_col
        mol_frame.mol_col = self.mol_col
        mol_frame.b64_col = self.b64_col


    def __getitem__(self, key):
        res = self.data[key]
        if isinstance(res, pd.DataFrame):
            result = self.new()
            result.data = res
            print_log(result.data, "subset")
        else:
            result = res
        return result


    def __setitem__(self, key, item):
        self.data[key] = item


    # def __getattr__(self, name):
    #     """Try to call undefined methods on the underlying pandas DataFrame."""
    #     def method(*args, **kwargs):
    #         res = getattr(self.data, name)(*args, **kwargs)
    #         if isinstance(res, pd.DataFrame):
    #             print("DataFrame")
    #             result = self.new()
    #             result.data = res
    #             print_log(result.data, name)
    #         else:
    #             result = res
    #         return result
    #     return method


    def print_log(self, component, add_info=""):
        if self.log:
            print_log(self.data, component, add_info)


    def __getattr__(self, name):
        """Try to call undefined methods on the underlying pandas DataFrame."""
        if hasattr(self.data, name):
            def method(*args, **kwargs):
                res = getattr(self.data, name)(*args, **kwargs)
                if isinstance(res, pd.DataFrame):
                    result = self.new()
                    result.data = res
                    print_log(result.data, name)
                else:
                    result = res
                return result
            return method
        else:
            raise AttributeError


    def _repr_html_(self):
        return self.data.to_html()


    def new(self):
        result = MolFrame()
        self._pass_properties(result)
        return result


    def copy(self):
        result = self.new()
        result.data = self.data.copy()
        return result


    def show(self, include_smiles=False, drop=[], **kwargs):
        if not include_smiles:
            drop.append(self.smiles_col)
        if len(drop) > 0:
            self.drop_cols(drop)
        return HTML(df_html(self.data))


    def view(self, title="MolFrame", include_smiles=False,
             drop=[], keep=[], fn="molframe.html", **kwargs):
        """Known kwargs: smiles_col, mol_col, id_col, selectable (bool), index (bool)"""
        self.add_mols()
        return view(self.data, title=title, include_smiles=include_smiles,
                    drop=drop, keep=keep, fn=fn,
                    smiles_col=self.smiles_col, mol_col=self.mol_col, _id_col=self.id_col,
                    **kwargs)


    def drop_cols(self, cols):
        """Drops the list of columns from the DataFrame.
        Listed columns that are not present in the DataFrame are simply ignored
        (no error is thrown)."""
        if self.inplace:
            drop_cols(self.data, cols, inplace=True)
            self.print_log("drop cols (inplace)")
        else:
            result = self.new()
            result.data = drop_cols(self.data, cols, inplace=False)
            print_log(result.data, "drop cols")
            return result


    def keep_cols(self, cols):
        if self.inplace:
            self.data = self.data[cols]
            self.print_log("keep cols (inplace)")
        else:
            result = self.new()
            result.data = self.data[cols]
            print_log(result.data, "keep cols")
            return result


    def load(self, fn, sep="\t"):
        """Read one or multiple result files and concatenate them into one MolFrame.
        `fn` is a single filename (string) or a list of filenames."""
        self.data = load(fn, sep=sep).data
        self.print_log("load data")


    def write_csv(self, fn, parameters=None, sep="\t"):
        result = self.data.copy()
        result = drop_cols(result, [self.mol_col])
        if isinstance(parameters, list):
            result = result[parameters]
        result.to_csv(fn, sep=sep, index=False)


    def write_pkl(self, fn):
        self.data.to_pickle(fn)


    def write_sdf(self, fn):
        self.add_mols()
        writer = Chem.SDWriter(fn)
        fields = [f for f in self.data.keys() if f != self.mol_col]
        for _, rec in self.data.iterrows():
            mol = rec[self.mol_col]
            for f in fields:
                if f in rec and rec[f]:
                    mol.SetProp(f, str(rec[f]))
            writer.write(mol)
        writer.close()


    def remove_mols(self):
        self.data.drop(self.mol_col, axis=1, inplace=True)


    def add_mols(self, force=False):
        # def _mol_from_smiles(smi):
        #     return pd.Series(mol_from_smiles(smi))
        if force or not self.has_mols:
            self.data[self.mol_col] = self.data[self.smiles_col].apply(mol_from_smiles)
            self.has_mols = True


    def add_smiles(self, isomeric_smiles=True):
        def _mol_to_smiles(mol):
            return pd.Series(Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles))
        self.data[self.smiles_col] = self.data[self.mol_col].apply(_mol_to_smiles)


    def b64_from_smiles(self):
        def _b64_from_smiles(smiles):
            mol = mol_from_smiles(smiles)
            result = b64.b64encode(pickle.dumps(mol)).decode()
            return result
        self.data[self.b64_col] = self.data[self.smiles_col].apply(_b64_from_smiles)


    def apply_to_mol(self, lambda_func, new_col_name):
        if self.inplace:
            self.data[new_col_name] = self.data[self.mol_col].apply(lambda_func)
        else:
            result = self.copy()
            result.data[new_col_name] = result.data[self.mol_col].apply(lambda_func)
            return result


    def mol_filter(self, query, add_h=False):
        query_mol = Chem.MolFromSmiles(query)
        if not query_mol:
            raise ValueError("Could not generate query mol.")
        if "[H]" in query or "#1" in query:
            add_h = True
            print("> explicit hydrogens turned on (add_h = True)")
        res_l = []
        for _, rec in self.data.iterrows():
            if self.mol_col in rec:
                mol = rec[self.mol_col]
            elif self.b64_col in rec:
                mol = pickle.loads(b64.b64decode(rec[self.b64_col]))
            elif self.smiles_col in rec:
                mol = Chem.MolFromSmiles(rec[self.smiles_col])
            else:
                raise ValueError("Neither Mol nor Mol_b64 nor Smiles column found.")
            if not mol: continue
            hit = False
            if add_h:
                mol_with_h = Chem.AddHs(mol)
                if mol_with_h.HasSubstructMatch(query_mol):
                    hit = True
            else:
                if mol.HasSubstructMatch(query_mol):
                    hit = True
            if hit:
                res_l.append(rec)
        result = self.new()
        result.data = pd.DataFrame(res_l)
        return result


def get_value(str_val):
    if not str_val:
        return None
    try:
        val = float(str_val)
        if "." not in str_val:
            val = int(val)
    except ValueError:
        val = str_val
    return val


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

    result.data = result.data.apply(pd.to_numeric, errors='ignore')
    print_log(result.data, "load molframe")
    return result


def load_sdf(fn):
    first_mol = True
    d = {"Mol": []}
    if isinstance(fn, str):
        file_obj = open(fn, "rb")
    else:
        file_obj = fn
    reader = Chem.ForwardSDMolSupplier(file_obj)
    for mol in reader:
        if first_mol:
            first_mol = False
            for prop in mol.GetPropNames():
                d[prop] = []
        # d["Smiles"].append(Chem.MolToSmiles(mol, isomericSmiles=True))
        for prop in mol.GetPropNames():
            if prop in d:
                d[prop].append(get_value(mol.GetProp(prop)))
            mol.ClearProp(prop)
        d["Mol"].append(mol)
    result = MolFrame()
    result.data = pd.DataFrame(d)
    result.has_mols = True
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


def load_resource(resource):
    res = resource.lower()
    glbls = globals()
    if "smi" in res:
        if "SMILES" not in glbls:
            # except NameError:
            print("- loading resource:                        (SMILES)")
            result = MolFrame()
            result.data = pd.read_csv(cprp.smiles_path, sep="\t")
            result.data = result.data[cprp.smiles_cols]
            result.data = result.data.apply(pd.to_numeric, errors='ignore')
            global SMILES
            SMILES = result
    elif "struct" in res:
        if "STRUCTURES" not in glbls:
            # except NameError:
            print("- loading resource:                        (STRUCTURES)")
            result = MolFrame()
            result.data = pd.read_csv(cprp.smiles_path, sep="\t")
            result.data = result.data[cprp.smiles_cols]
            result.data = result.data.apply(pd.to_numeric, errors='ignore')
            result.data[result.mol_col] = result.data["Mol_b64"].apply(pickle.loads(b64.b64decode))
            result.data.drop("Mol_b64", axis=1, inplace=True)
            global STRUCT
            STRUCT = result
    elif "cont" in res:
        if "CONTAINER" not in glbls:
            print("- loading resource:                        (CONTAINER)")
            result = pd.read_csv(cprp.container_data_path, sep="\t")
            if len(cprp.container_data_cols) > 0:
                result = result[cprp.container_data_cols]
            result = result.apply(pd.to_numeric, errors='ignore')
            global CONTAINER
            CONTAINER = MolFrame()
            CONTAINER.data = result
    elif "batch" in res:
        if "BATCH" not in glbls:
            print("- loading resource:                        (BATCH)")
            result = pd.read_csv(cprp.batch_data_path, sep="\t")
            if len(cprp.batch_data_cols) > 0:
                result[cprp.batch_data_cols]
            result = result.apply(pd.to_numeric, errors='ignore')
            global BATCH
            BATCH = MolFrame()
            BATCH.data = result
    else:
        raise FileNotFoundError("# unknown resource: {}".format(resource))
