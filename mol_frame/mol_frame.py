#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
########
MolFrame
########

*Created on Thu Aug 31, 2017 by A. Pahl*

Tools for working with MolFrame dataframes."""

# import pdb

import base64 as b64
import gzip
import os.path as op
import pickle  # for b64 encoded mol objects
import sys
import tempfile
import time
from copy import deepcopy

import numpy as np
import pandas as pd
import rdkit.Chem.Descriptors as Desc
from rdkit import DataStructs
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.MolStandardize.charge import Uncharger
from rdkit.Chem.MolStandardize.fragment import LargestFragmentChooser
from rdkit.Chem.MolStandardize.standardize import Standardizer

from mol_frame import nb_tools as nbt
from mol_frame import templ
from mol_frame import tools as mft
from mol_frame.mol_images import b64_mol, check_2d_coords
from mol_frame.viewers import mol_img_tag, show_grid, write_grid

molvs_s = Standardizer()
molvs_l = LargestFragmentChooser()
molvs_u = Uncharger()



x = b64_mol  # make the linter shut up

try:
    from misc_tools import apl_tools

    AP_TOOLS = True
    #: Library version
    VERSION = apl_tools.get_commit(__file__)
    # I use this to keep track of the library versions I use in my project notebooks
    print("{:45s} ({})".format(__name__, VERSION))

except ImportError:
    AP_TOOLS = False
    print(
        "{:45s} ({})".format(
            __name__,
            time.strftime("%y%m%d-%H:%M", time.localtime(op.getmtime(__file__))),
        )
    )

try:
    import holoviews as hv

    hv.extension("bokeh")
    from bokeh.models import HoverTool

    HOLOVIEWS = True

except ImportError:
    HOLOVIEWS = False
    print(
        "* holoviews could not be imported. scatter() and struct_hover() are not available."
    )


IPYTHON = nbt.is_interactive_ipython()

if IPYTHON:
    from IPython.core.display import HTML


DEBUG = False
nbits = 1024
FPDICT = {}
FPDICT["ecfp0"] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 0, nBits=nbits)
FPDICT["ecfp2"] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 1, nBits=nbits)
FPDICT["ecfp4"] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 2, nBits=nbits)
FPDICT["ecfp6"] = lambda m: Chem.GetMorganFingerprintAsBitVect(m, 3, nBits=nbits)
FPDICT["ecfc0"] = lambda m: Chem.GetMorganFingerprint(m, 0)
FPDICT["ecfc2"] = lambda m: Chem.GetMorganFingerprint(m, 1)
FPDICT["ecfc4"] = lambda m: Chem.GetMorganFingerprint(m, 2)
FPDICT["ecfc6"] = lambda m: Chem.GetMorganFingerprint(m, 3)
FPDICT["fcfp2"] = lambda m: Chem.GetMorganFingerprintAsBitVect(
    m, 1, useFeatures=True, nBits=nbits
)
FPDICT["fcfp4"] = lambda m: Chem.GetMorganFingerprintAsBitVect(
    m, 2, useFeatures=True, nBits=nbits
)
FPDICT["fcfp6"] = lambda m: Chem.GetMorganFingerprintAsBitVect(
    m, 3, useFeatures=True, nBits=nbits
)
FPDICT["fcfc2"] = lambda m: Chem.GetMorganFingerprint(m, 1, useFeatures=True)
FPDICT["fcfc4"] = lambda m: Chem.GetMorganFingerprint(m, 2, useFeatures=True)
FPDICT["fcfc6"] = lambda m: Chem.GetMorganFingerprint(m, 3, useFeatures=True)


def debug_print(txt, val):
    if DEBUG:
        txt = txt + ":"
        print("DEBUG   {:20s}".format(txt), val)


class MolFrame(object):
    """A wrapper class around Pandas DataFrame with added structure handling capabilities."""

    def __init__(self, init=None, log=True):
        if init is None:
            self.data = None
        else:
            self.data = pd.DataFrame(init)
        self.inplace = False
        self.has_mols = False
        self.id_col = "Compound_Id"
        self.smiles_col = "Smiles"
        self.mol_col = "Mol"
        self.b64_col = "Mol_b64"
        self.fp_col = "FP_b64"
        self.fp_name = ""  # Fingerprint method, e.g. ecfc4

    def _pass_properties(self, mol_frame):
        mol_frame.inplace = self.inplace
        mol_frame.has_mols = self.has_mols
        mol_frame.id_col = self.id_col
        mol_frame.smiles_col = self.smiles_col
        mol_frame.mol_col = self.mol_col
        mol_frame.b64_col = self.b64_col
        mol_frame.fp_col = self.fp_col
        mol_frame.fp_name = self.fp_name

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
        drop = [self.b64_col, self.mol_col, "Image"]
        result = drop_cols(self.data, drop)
        return result.to_html()

    def new(self):
        result = MolFrame()
        self._pass_properties(result)
        return result

    def copy(self):
        result = self.new()
        result.data = self.data.copy()
        return result

    def mol_img_tag(self, x):
        if self.use_col is None:
            self.find_mol_col()
        mol = self.mol_method(x)
        return mol_img_tag(mol)

    def to_html(self, **kwargs):
        """New table display, relying on pd.to_html and formatters.
        The `Molecule` column is put first.

        Known kwargs:
            formatters (dict or None)
            escape (bool): Default is False.
            index (bool): Whether to display the index or not, default is False.
            rename_mol_col (bool): Default is True.
            selectable (bool): Default is False."""
        self.find_mol_col()
        formatters = {self.use_col: self.mol_img_tag}
        formatters.update(kwargs.pop("formatters", {}))
        escape = kwargs.pop("escape", False)
        index = kwargs.pop("index", False)
        rename_mol_col = kwargs.pop("rename_mol_col", True)
        selectable = kwargs.pop("selectable", False)
        columns = list(self.data.columns)
        columns = [columns.pop(columns.index(self.use_col))] + columns
        tmp = self.data.copy()
        if selectable:
            columns = ["$Sel$"] + columns
            tmp["$Sel$"] = ""
        tmp = tmp[columns]
        result = tmp.to_html(formatters=formatters, escape=escape, index=index, **kwargs)

        if rename_mol_col:
            result = result.replace(f">{self.use_col}</th>", ">Molecule</th>")
        return result

    def show(self, **kwargs):
        """New table display, relying on pd.to_html and formatters.

        Known kwargs:
            formatters (dict or None)
            escape (bool): Default is False.
            index (bool): Whether to display the index or not, default is False.
            rename_mol_col (bool)"""
        kwargs.pop("selectable", False)
        return HTML(self.to_html(**kwargs))

    def show_grid(
        self,
        include_smiles=False,
        title="MolFrame",
        drop=[],
        keep=[],
        smiles_col="Smiles",
        mol_col="Mol",
        id_col=None,
        b64_col=None,
        fp_col=None,
        **kwargs
    ):
        self.add_mols()
        return show_grid(
            self.data,
            title=title,
            include_smiles=include_smiles,
            drop=drop,
            keep=keep,
            smiles_col=self.smiles_col,
            mol_col=self.mol_col,
            id_col=self.id_col,
            b64_col=self.b64_col,
            fp_col=self.fp_col,
            **kwargs
        )

    def write_tbl(
        self,
        title="MolFrame",
        fn="molframe.html",
        **kwargs
    ):
        """Known kwargs:
            selectable (bool); Currently without function
            index (bool): Whether to show the index or not, default is False.
            intro (str): Optional HTML text that will be inserted above the table.
            format (str): Formatting used for the table. Available options: ``simple``, ``bootstrap``.
            formatters (dict or None)
            escape (bool): Default is False.
            rename_mol_col (bool): Default is True.
            selectable (bool): Default is False."""
        tbl_format = kwargs.pop("format", "bootstrap").lower()
        intro = kwargs.pop("intro", "")
        selectable = kwargs.get("selectable", False)
        if "boot" in tbl_format:
            # Bootstrap does not seem to play nicely with the pandas index:
            kwargs["index"] = False
        table = self.to_html(**kwargs)
        d = {}
        if "boot" in tbl_format:
            table = templ.bootstrap_options(table, selectable, id_col=self.id_col)
            if selectable:
                t = templ.MFTemplate(templ.TABLE_BOOTSTRAP_SELECT)
                d["id_col"] = self.id_col + "s"  # Mehrzahl
            else:
                t = templ.MFTemplate(templ.TABLE_BOOTSTRAP)
        else:
            print("* simple format used.")
            t = templ.MFTemplate(templ.TABLE_SIMPLE)
        d.update({"title": title, "intro": intro, "table": table})
        page = t.substitute(d)
        templ.write(page, fn)
        if IPYTHON:
            return HTML(f'<a href="{fn}">{title}</a>')

    def write_grid(
        self, title="MolGrid", fn="molgrid.html", **kwargs
    ):
        """Known kwargs: interactive (bool)
                         highlight (dict)
                         header (str)
                         summary (str)
                         hlsss (colname)
                         interactive (bool)
                         link_templ. link_col (str) (then interactive is false)

        Example:
            title = f"Sim Scaf {scaf_num}"
            fn = f"chembl_sim_scaf_{scaf_num}.html"
            header = f"ChEMBL Compounds Most Similar to Scaf {scaf_num}"
            df_sim = mf.read_csv(f"chembl_sim_scaf_{scaf_num}.tsv")
            df_sim.id_col = "Chembl_Id"
            df_sim = df_sim.sort_values("Sim", ascending=False)
            df_sim = df_sim.head(100)
            img_tag = mfi.mol_img_tag(smiles[scaf_num - 1])
            highest_sim = df_sim["Sim"].values[0]
            summary = f'''Compounds from ChEMBL 24 with similarity to scaffold {scaf_num}. <br>{img_tag}<br><br>
            The highest observed similarity was {highest_sim} %.<br>
            The search was performed using Morgan fingerprints and Tanimoto similarity
            on the Murcko scaffolds of the compounds from ChEMBL 24).<br>'''
            link_col = "Chembl_Id"
            link_templ = "https://www.ebi.ac.uk/chembl/compound/inspect/{}"
            """
        self.add_mols()
        if self.id_col is not None and self.id_col not in self.data.keys():
            self.id_col = None
        return write_grid(
            self.data,
            title=title,
            fn=fn,
            smiles_col=self.smiles_col,
            mol_col=self.mol_col,
            id_col=self.id_col,
            b64_col=self.b64_col,
            fp_col=self.fp_col,
            **kwargs
        )

    def grid(self, title="MolGrid", drop=[], keep=[], fn="molgrid.html", **kwargs):
        """Known kwargs: interactive (bool)
                         highlight (dict)
                         hlsss (colname)
                         truncate (int)"""
        self.add_mols()
        if self.id_col is not None and self.id_col not in self.data.keys():
            self.id_col = None
        return show_grid(
            self.data,
            title=title,
            drop=drop,
            keep=keep,
            fn=fn,
            smiles_col=self.smiles_col,
            mol_col=self.mol_col,
            id_col=self.id_col,
            b64_col=self.b64_col,
            fp_col=self.fp_col,
            **kwargs
        )

    def info(self):
        """Show a summary of the MolFrame."""
        keys = list(self.data.columns)
        info = []
        for k in keys:
            info.append(
                {
                    "Field": k,
                    "Count": self.data[k].notna().count(),
                    "Type": str(self.data[k].dtype),
                }
            )
        info.append({"Field": "Total", "Type": "", "Count": self.data.shape[0]})
        return pd.DataFrame(info)

    def compute(self):
        """Kept for backwards compatibility."""
        df = self.data.copy()
        result = self.new()
        result.data = df
        print_log(df, "compute")
        return result

    def groupby(self, by=None, num_agg=["median", "mad", "count"], str_agg="unique"):
        if by is None:
            by = self.id_col
        result = self.new()
        result.data = groupby(self.data, by=by, num_agg=num_agg, str_agg=str_agg)
        print_log(self.data, "groupby")
        return result

    def concat(self, other):
        if hasattr(other, "data"):  # a MolFrame instance
            df = other.data
        else:
            df = other
        result = self.new()
        result.data = pd.concat([self.data, df])
        print_log(result.data, "concat")
        return result

    def keep_cols(self, cols):
        """Keeps the list of columns of the DataFrame.
        Listed columns that are not present in the DataFrame are simply ignored
        (no error is thrown)."""
        if self.inplace:
            keep_cols(self.data, cols, inplace=True)
            self.print_log("keep cols (inplace)")
        else:
            result = self.new()
            result.data = keep_cols(self.data, cols, inplace=False)
            print_log(result.data, "keep cols")
            return result

    def drop_cols(self, cols):
        """Drops the list of columns from the DataFrame.
        Listed columns that are not present in the DataFrame are simply ignored
        (no error is thrown)."""
        if isinstance(cols, str):
            cols = [cols]
        if self.mol_col in cols:
            self.has_mols = False
        if self.inplace:
            drop_cols(self.data, cols, inplace=True)
            self.print_log("drop cols (inplace)")
        else:
            result = self.new()
            result.data = drop_cols(self.data, cols, inplace=False)
            print_log(result.data, "drop cols")
            return result

    def read_csv(self, fn, sep="\t"):
        """Read one or multiple result files and concatenate them into one MolFrame.
        ``fn`` is a single filename (string) or a list of filenames."""
        self.data = read_csv(fn, sep=sep).data
        self.print_log("read CSV")

    def write_csv(self, fn, parameters=None, sep="\t"):
        result = self.data
        result = drop_cols(result, [self.mol_col])
        if isinstance(parameters, list):
            result = result[parameters]
        result.to_csv(fn, sep=sep, index=False)

    def write_tmp(self, parameters=None):
        """Write the Mol_Frame to a temporary file. Returns the filename.
        Can be re-openend with mf.read_csv(fn)."""
        if self.is_dask:  # Dask
            glob_str = "-*"
        else:
            glob_str = ""  # single file
        fn = op.join(
            tempfile.gettempdir(),
            "mf_{}{}.tmp".format(time.strftime("%Y-%m-%d_%H%M%S"), glob_str),
        )
        self.write_csv(fn, parameters=parameters, sep="\t")
        return fn

    def write_pkl(self, fn):
        """Only works when the data object is a Pandas DataFrame."""
        pkl = [
            self.inplace,
            self.has_mols,
            self.id_col,
            self.smiles_col,
            self.mol_col,
            self.b64_col,
            self.fp_col,
            self.fp_name,
            self.data,
        ]
        with open(fn, "wb") as f:
            pickle.dump(pkl, f)

    def write_sdf(self, fn):
        """Only works when the data object is a Pandas DataFrame."""

        writer = Chem.SDWriter(fn)
        fields = []
        for f in self.data.keys():
            if f != self.mol_col and f != self.b64_col:
                fields.append(f)
        self.find_mol_col()
        for _, rec in self.data.iterrows():
            mol = self.mol_method(rec[self.use_col])
            for f in fields:
                if f in rec and rec[f]:
                    mol.SetProp(f, str(rec[f]))
            writer.write(mol)
        writer.close()

    def remove_mols(self):
        """Remove the mol objects from the MolFrame."""
        if self.inplace:
            self.data.drop(self.mol_col, axis=1, inplace=True)
            self.has_mols = False
        else:
            result = self.new()
            result.data = self.data.drop(self.mol_col, axis=1)
            result.has_mols = False
            print_log(result.data, "remove mols")
            return result

    def remove_smiles_and_b64(self):
        """Remove the Smiles and Mol_b64 columns from the MolFrame, if present."""
        drop = []
        for k in [self.smiles_col, self.b64_col]:
            if k in self.data.columns:
                drop.append(k)
        if self.inplace:
            self.data.drop(drop, axis=1, inplace=True)
        else:
            result = self.new()
            result.data = self.data.drop(drop, axis=1)
            print_log(result.data, "drop cols")
            return result

    def add_mols(self, force=False, standardize=True, remove_src=False):
        """Adds mol objects to the MolFrame.
        Tries to use the Mol_b64 column. If that is not present, it uses Smiles.
        Only works when the data object is a Pandas DataFrame.
        Standardizes by default (see standardize_mols())."""
        self.find_mol_col()
        if self.inplace:
            if force or not self.has_mols:
                self.data[self.mol_col] = self.data[self.use_col].apply(self.mol_method)
                self.has_mols = True
                if remove_src:
                    self.data.drop(self.use_col, axis=1, inplace=True)
                if standardize:
                    self.standardize_mols()
                print_log(self.data, "add mols")
        else:
            result = self.copy()
            if force or not self.has_mols:
                result.data[self.mol_col] = result.data[self.use_col].apply(
                    self.mol_method
                )
                result.has_mols = True
                if remove_src:
                    result.data = result.data.drop(self.use_col, axis=1)
            if standardize:
                result = result.standardize_mols()
                print_log(result.data, "add mols")
            return result

    def add_images(self, force=False):
        """Adds an Image column to the MolFrame, used for structure tooltips in plotting.
        Only works on Pandas DataFrames, does not work for Dask DataFrames
        (call `.compute()` first)."""
        if "Image" in self.keys() and not force:
            if self.inplace:
                return
            else:
                result = self.new()
                result.data = self.data
                print_log(result.data, "add img")
                return result

        self.find_mol_col()

        def _img_method(x):
            return "data:image/png;base64,{}".format(b64_mol(self.mol_method(x)))

        if self.inplace:
            self.data["Image"] = self.data[self.use_col].apply(_img_method)
        else:
            result = self.new()
            result.data = self.data
            result.data["Image"] = result.data[self.use_col].apply(_img_method)
            print_log(result.data, "add img")
            return result

    def add_smiles(self, isomeric_smiles=True):
        """Adds Smiles column from mol objects to the MolFrame."""

        def _smiles_from_mol(mol):
            return Chem.MolToSmiles(mol, isomericSmiles=isomeric_smiles)

        if self.inplace:
            self.apply_to_mol(self.smiles_col, _smiles_from_mol)
        else:
            result = self.apply_to_mol(self.smiles_col, _smiles_from_mol)
            print_log(result.data, "add smiles")
            return result

    def standardize_mols(self):
        """Standardizes the molecules in the Mol column.
        Applies MolVS Standardizer, LargestFragment and Uncharger.
        Requires the MolFrame to have mols (has_mols == True)."""
        if not self.has_mols:
            print("* MolFrame does not have mols. Run add_mols() first.")
            return

        if self.inplace:
            self.data[self.mol_col] = self.data[self.mol_col].apply(standardize_mol)
        else:
            result = self.copy()
            result.data[self.mol_col] = result.data[self.mol_col].apply(standardize_mol)
            print_log(result.data, "add smiles")
            return result

    def add_inchikeys(self):
        """Adds Inchi Keys."""
        self.find_mol_col()
        if len(self.data) > 5000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False

        def _lambda(x):
            if show_prog:
                pb.inc()
            mol = self.mol_method(x)
            if not mol:
                return "NO_MOL."
            try:
                ik = Chem.inchi.MolToInchiKey(mol)
            except ValueError:
                ik = "FAILED."
            return ik

        if self.inplace:
            self.data["InchiKey"] = self.data[self.use_col].apply(_lambda)
            if show_prog:
                pb.done()
        else:
            result = self.copy()
            result.data["InchiKey"] = result.data[self.use_col].apply(_lambda)
            if show_prog:
                pb.done()
            return result

    def add_fp(self, fp_name="ecfc4"):
        """Adds a fingerprint column to the MolFrame.
        The available fingerprints can be viewed by:
        >>> print(mf.FDICT.keys())"""

        fp_name = fp_name.lower()

        def fp_method(x):
            return b64.b64encode(pickle.dumps(FPDICT[fp_name](x))).decode()

        if self.inplace:
            self.fp_name = fp_name
            self.apply_to_mol(self.fp_col, fp_method)
        else:
            result = self.apply_to_mol(self.fp_col, fp_method)
            result.fp_name = fp_name
            print_log(result.data, "add fp")
            return result

    def add_b64(self):
        """Adds Mol_b64 column to MolFrame."""

        def _b64_from_mol(mol):
            result = b64.b64encode(pickle.dumps(mol)).decode()
            return result

        if self.inplace:
            self.apply_to_mol(self.b64_col, _b64_from_mol)
        else:
            result = self.apply_to_mol(self.b64_col, _b64_from_mol)
            print_log(result.data, "add b64")
            return result

    def find_mol_col(self):
        """Find a suitable mol column.
        Either Mol, Mol_b64 or Smiles."""
        self.use_col = None
        if self.has_mols:
            self.use_col = self.mol_col
            self.mol_method = lambda x: x
        elif self.b64_col in self.data.columns:
            print("* using", self.b64_col)
            self.use_col = self.b64_col
            self.mol_method = lambda x: pickle.loads(b64.b64decode(x))
        elif self.smiles_col in self.data.columns:
            print("* using", self.smiles_col)
            self.use_col = self.smiles_col
            self.mol_method = lambda x: mol_from_smiles(x)
        else:
            raise KeyError("No suitable Mol column found.")

    def apply_to_col(self, col_name, new_col_name, lambda_func):
        """Applies a func to a column in the MolFrame.
        A wrapper around pd.apply to enable progress bars.
        Returns a new copy or modifies inplace, depending on self.inplace."""
        if len(self.data) > 5000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False

        def _apply(x):
            if show_prog:
                pb.inc()
            return lambda_func(x)

        if self.inplace:
            self.data[new_col_name] = self.data[col_name].apply(_apply)
            if show_prog:
                pb.done()
        else:
            result = self.new()
            result.data = self.data
            result.data[new_col_name] = result.data[col_name].apply(_apply)
            if show_prog:
                pb.done()
            return result

    def apply_to_mol(self, new_col_name, lambda_func):
        """Applies a func to the Mol object, which is generated on-the-fly, if necessary.
        Displays a progress bar for longer operations.
        Returns a new copy or modifies inplace, depending on self.inplace."""
        self.find_mol_col()
        if len(self.data) > 1000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False

        def _apply(x):
            if show_prog:
                pb.inc()
            mol = self.mol_method(x)
            if not mol:
                return pd.np.nan
            return lambda_func(mol)

        if self.inplace:
            self.data[new_col_name] = self.data[self.use_col].apply(_apply)
            if show_prog:
                pb.done()
        else:
            result = self.new()
            result.data = self.data.copy()
            result.data[new_col_name] = result.data[self.use_col].apply(_apply)
            if show_prog:
                pb.done()
            return result

    def apply_numeric(self):
        """Apply pd.numeric."""
        if self.inplace:
            self.data.apply(pd.to_numeric, errors="ignore", axis=1)
        else:
            result = self.new()
            result.data = self.data.apply(pd.to_numeric, errors="ignore", axis=1)
            return result

    def check_2d_coords(self, force=False):
        """Generates 2D coordinates if necessary.
        Requires the Mol object to be present (use add_mols() )."""
        self.find_mol_col()
        if len(self.data) > 1000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False

        def _apply(x):
            if show_prog:
                pb.inc()
            mol = self.mol_method(x)
            if mol:
                check_2d_coords(mol, force=force)

        if self.inplace:
            self.data[self.use_col].apply(_apply)
            if show_prog:
                pb.done()
        else:
            result = self.copy()
            result.use_col = self.use_col
            result.mol_method = self.mol_method
            result.data[self.use_col].apply(_apply)
            if show_prog:
                pb.done()
            return result

    def rescale(self, f=1.5):
        def _transform(m):
            if show_prog:
                pb.inc()
            tm = np.zeros((4, 4), np.double)
            for i in range(3):
                tm[i, i] = f
            tm[3, 3] = 1.0
            Chem.TransformMol(m, tm)

        self.find_mol_col()
        if len(self.data) > 1000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False

        if self.inplace:
            if not self.has_mols:
                return
            self.data[self.use_col].apply(_transform)
        else:
            result = self.copy()
            result.use_col = self.use_col
            result.mol_method = self.mol_method
            if not self.has_mols:
                return result
            result.data[self.use_col].apply(_transform)
            if show_prog:
                pb.done()
            return result

    def keep_largest_fragment(self):
        """Removes salts, etc.
        Returns the new molecules as SmilesFrag Column."""

        def _largest_frag(mol):
            mols = Chem.GetMolFrags(mol, asMols=True)
            if len(mols) > 1:
                frag_ctr[0] += 1
                mols = sorted(mols, key=Desc.HeavyAtomCount, reverse=True)
                new_smiles = Chem.MolToSmiles(mols[0], isomericSmiles=True)
            else:
                new_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            return new_smiles

        frag_ctr = [0]
        if self.inplace:
            self.apply_to_mol("SmilesFrag", _largest_frag)
            print("* fragments removed in {} molecules.".format(frag_ctr[0]))
        else:
            result = self.apply_to_mol("SmilesFrag", _largest_frag)
            print("* fragments removed in {} molecules.".format(frag_ctr[0]))
            return result

    def id_filter(self, cpd_ids, reset_index=True, sort_by_input=False):
        """Returns a new MolFrame instance consisting only of the given cpd_ids.
        When ``sort_by_input == True``, the lines in the new MolFrame will have the same order
        as the given cpd_ids."""
        if not isinstance(cpd_ids, list):
            cpd_ids = [cpd_ids]
        df = self.data[self.data[self.id_col].isin(cpd_ids)]

        if reset_index:
            df.reset_index(inplace=True)
            df.drop("index", axis=1, inplace=True)
        if sort_by_input:  # requires the data obejct to be a Pandas DataFrame.
            if self.is_dask:
                print("  - computing Pandas DataFrame...")
                df = df.compute()
            df["_sort"] = pd.Categorical(
                df[self.id_col], categories=cpd_ids, ordered=True
            )
            df = df.sort_values("_sort")
            df.drop("_sort", axis=1, inplace=False)
        result = self.new()
        result.data = df
        print_log(result.data, "id filter")
        return result

    def mol_filter(self, query, add_h=False):
        """Substructure filter. Returns a new MolFrame instance.
        ``query`` has to be a Smiles string."""
        if len(self.data) > 5000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False
        query_mol = Chem.MolFromSmiles(query)
        if not query_mol:
            raise ValueError("Could not generate query mol.")
        if "[H]" in query or "#1" in query:
            add_h = True
            print("> explicit hydrogens turned on (add_h = True)")
        res_l = []
        self.find_mol_col()
        for _, rec in self.data.iterrows():
            if show_prog:
                pb.inc()
            mol = self.mol_method(rec[self.use_col])
            if not mol:
                continue
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
        if show_prog:
            pb.done()
        print_log(result.data, "mol_filter")
        return result

    def sim_filter(self, query, cutoff=0.75):
        """Similarity filter. Returns a new MolFrame instance.
        Add a suitable fingerprint once with addf_fps(),
        then give a reference molecule or a SMILES string as query."""
        if len(self.fp_name) == 0 or self.fp_col not in self.data.columns:
            raise KeyError(
                "No fingerprints found. Please generate them first with add_fp()."
            )
        if len(self.data) > 5000:
            show_prog = True
            pb = nbt.Progressbar(end=len(self.data))
        else:
            show_prog = False
        if isinstance(query, str):
            query_mol = Chem.MolFromSmiles(query)
        else:
            query_mol = deepcopy(query)
        if not query_mol:
            raise ValueError("Could not generate query mol.")
        fp_method = FPDICT[self.fp_name]
        query_fp = fp_method(query_mol)
        res_l = []
        for _, rec in self.data.iterrows():
            if show_prog:
                pb.inc()
            mol_fp = pickle.loads(b64.b64decode(rec[self.fp_col]))
            sim = DataStructs.TanimotoSimilarity(query_fp, mol_fp)
            if sim >= cutoff:
                rec["Sim"] = sim
                res_l.append(rec)
        result = self.new()
        result.data = pd.DataFrame(res_l)
        print_log(result.data, "sim_filter")
        if show_prog:
            pb.done()
        return result

    def scatter(
        self,
        x,
        y,
        colorby=None,
        options={},
        styles={},
        title="Scatter Plot",
        force=False,
    ):
        """Possible options: width, height, legend_position [e.g. "top_right"]
        Possible styles: size, cmap [brg, Accent, rainbow, jet, flag, Wistia]
        Only works when the data object is a Pandas DataFrame."""
        if not HOLOVIEWS:
            print("* HoloViews not available.")
            return None
        hover = struct_hover(self, force=force)
        plot_options = {
            "width": 800,
            "height": 450,
            "legend_position": "top_left",
            "tools": [hover],
        }
        plot_styles = {"size": 8, "cmap": "brg"}
        plot_options.update(options)
        plot_styles.update(styles)
        vdims = [y, self.id_col, "Image"]
        if colorby is not None:
            vdims.append(colorby)
            plot_options["color_index"] = len(vdims)
        opts = {"Scatter": {"plot": plot_options, "style": plot_styles}}
        scatter_plot = hv.Scatter(self.data, x, vdims=vdims, label=title)
        return scatter_plot(opts)


def groupby(df_in, by=None, num_agg=["median", "mad", "count"], str_agg="unique"):
    """Other str_aggs: "first", "unique"."""

    def _concat(values):
        return "; ".join(str(x) for x in values)

    def _unique(values):
        return "; ".join(set(str(x) for x in values))

    if isinstance(num_agg, str):
        num_agg = [num_agg]
    df_keys = df_in.columns
    numeric_cols = list(df_in.select_dtypes(include=[np.number]).columns)
    str_cols = list(set(df_keys) - set(numeric_cols))
    # if by in numeric_cols:
    try:
        by_pos = numeric_cols.index(by)
        numeric_cols.pop(by_pos)
    except ValueError:
        pass
    try:
        by_pos = str_cols.index(by)
        str_cols.pop(by_pos)
    except ValueError:
        pass
    aggregation = {}
    for k in numeric_cols:
        aggregation[k] = num_agg
    if str_agg == "join":
        str_agg_method = _concat
    elif str_agg == "first":
        str_agg_method = "first"
    elif str_agg == "unique":
        str_agg_method = _unique
    for k in str_cols:
        aggregation[k] = str_agg_method
    df = df_in.groupby(by)
    df = df.agg(aggregation).reset_index()
    df_cols = [
        "_".join(col).strip("_").replace("_<lambda>", "").replace("__unique", "")
        for col in df.columns.values
    ]
    df.columns = df_cols
    return df


def get_value(str_val):
    """convert a string into float or int, if possible."""
    if not str_val:
        return None
    try:
        val = float(str_val)
        if "." not in str_val:
            val = int(val)
    except ValueError:
        val = str_val
    return val


def print_log(df, component, add_info=""):
    component = component + ":"
    if len(add_info) > 0:
        add_info = "    ({})".format(add_info)
    if hasattr(df, "shape"):
        print(
            "* {:22s} ({:5d} | {:4d}){}".format(
                component, df.shape[0], df.shape[1], add_info
            )
        )
    else:
        print("* {:22s} (Dask  | {:4d}){}".format(component, len(df.columns), add_info))


def read_csv(fn, sep="\t"):
    """Read one or multiple result files and concatenate them into one MolFrame.
    ``fn`` is a single filename (string) or a list of filenames."""
    result = MolFrame()
    if isinstance(fn, list):
        df_list = []
        for f in fn:
            df_list.append(pd.read_csv(f, sep=sep))
        result.data = pd.concat(df_list)
    else:
        result.data = pd.read_csv(fn, sep=sep)
    # result.data = result.data.apply(pd.to_numeric, errors='ignore', axis=1)
    print_log(result.data, "read CSV")
    return result


def read_sdf(fn, store_mol_as="Mol_b64", gen2d=False):
    """Create a MolFrame instance from an SD file (can be gzipped (fn ends with ``.gz``)).

    Arguments:
        store_mol_as: "Mol_b64" or "Smiles" """
    if store_mol_as not in ["Mol_b64", "Smiles"]:
        print("* Mols are stored as Mol_b64")
        store_mol_as = "Mol_b64"
    d = {}
    ctr = 0
    first_mol = True
    if not isinstance(fn, list):
        fn = [fn]
    for f in fn:
        do_close = True
        if isinstance(f, str):
            if f.endswith(".gz"):
                file_obj = gzip.open(f, mode="rb")
            else:
                file_obj = open(f, "rb")
        else:
            file_obj = f
            do_close = False
        reader = Chem.ForwardSDMolSupplier(file_obj)
        for mol in reader:
            ctr += 1
            if not mol:
                print(ctr, file=sys.stderr)
                continue
            if first_mol:
                first_mol = False
                for prop in mol.GetPropNames():
                    if prop in [store_mol_as, "order"]:
                        continue
                    d[prop] = []
                d_keys = set(d.keys())
                d[store_mol_as] = []
            mol_props = set()
            for prop in mol.GetPropNames():
                if prop in d_keys:
                    mol_props.add(prop)
                    d[prop].append(get_value(mol.GetProp(prop)))
                mol.ClearProp(prop)

            # append NAN to the missing props that were not in the mol:
            missing_props = d_keys - mol_props
            for prop in missing_props:
                d[prop].append(np.nan)
            if gen2d:
                check_2d_coords(mol, force=True)
            if store_mol_as == "Smiles":
                smi = Chem.MolToSmiles(mol)
                d["Smiles"].append(smi)
            else:
                mol_b64 = b64.b64encode(pickle.dumps(mol)).decode()
                d["Mol_b64"].append(mol_b64)
        if do_close:
            file_obj.close()
    for k in d.keys():
        print(len(d[k]))
    result = MolFrame()
    result.data = pd.DataFrame(d)
    print_log(result.data, "read SDF")
    return result


def read_pkl(fn):
    with open(fn, "rb") as f:
        pkl = pickle.load(f)
    result = MolFrame()
    result.inplace = pkl[0]
    result.has_mols = pkl[1]
    result.id_col = pkl[2]
    result.smiles_col = pkl[3]
    result.mol_col = pkl[4]
    result.b64_col = pkl[5]
    result.fp_col = pkl[6]
    result.fp_name = pkl[7]
    result.data = pkl[8]
    return result


def mol_from_smiles(smi, calc_2d=True):
    """Generate a mol from Smiles.
    For invalid Smiles it generates a No-structure."""
    if smi == "foo":  # dask artifact
        smi = "*"
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles("*")
    else:
        if calc_2d:
            check_2d_coords(mol)
    return mol


def keep_cols(df, cols, inplace=False):
    """Keep the list of columns of the DataFrame.
    Listed columns that are not present in the DataFrame are simply ignored
    (no error is thrown)."""
    if isinstance(cols, str):
        cols = [cols]
    if inplace:
        drop = list(set(df.columns) - set(cols))
        df.drop(drop, axis=1, inplace=True)
    else:
        keep = list(set(cols).intersection(df.columns))
        result = df[keep]
        return result


def drop_cols(df, cols, inplace=False):
    """Drops the list of columns from the DataFrame.
    Listed columns that are not present in the DataFrame are simply ignored
    (no error is thrown)."""
    if isinstance(cols, str):
        cols = [cols]
    drop = set(cols).intersection(set(df.columns))
    if inplace:
        df.drop(drop, axis=1, inplace=True)
    else:
        result = df.drop(drop, axis=1)
        return result


def load_resources():
    """Known resources:
        SMILES,
        STRUCTURES: containing Mol_b64 column,
        BATCH, CONTAINER, DATA"""
    mf_config = mft.load_config("config")
    global SMILES
    SMILES = mf_config["Paths"]["SmilesPath"]
    global STRUCT
    STRUCT = mf_config["Paths"]["StructPath"]
    global DATA
    DATA = mf_config["Paths"]["ContainerDataPath"]
    global CONTAINER
    CONTAINER = mf_config["Paths"]["ContainerPath"]
    global BATCH
    BATCH = mf_config["Paths"]["BatchPath"]


def struct_hover(mf, force=False):
    """Create a structure tooltip that can be used in Holoviews.
    Takes a MolFrame instance as parameter."""
    if not HOLOVIEWS:
        print("* HoloViews not available.")
        return None
    mf.add_images(force=force)
    hover = HoverTool(
        tooltips="""
            <div>
                <div>
                    <img src="@Image" alt="Mol" width="70%"><br>
                <div>
                <div>
                    <span style="font-size: 12px; font-weight: bold;">@{}</span>
                </div>
            </div>
        """.format(
            mf.id_col
        )
    )
    return hover


def standardize_mol(mol):
    if mol is None: return None
    mol = molvs_s.standardize(mol)
    mol = molvs_l.choose(mol)
    mol = molvs_u.uncharge(mol)
    return mol
