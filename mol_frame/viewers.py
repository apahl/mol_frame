#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###########################
MolFrame Viewer Components
###########################

*Created on Wed Aug 30 9:00 2017 by A. Pahl*

View MolFrames in the Browser.
Code is based on bluenote10's NimData html browser
(https://github.com/bluenote10/NimData/blob/master/src/nimdata/html.nim)."""

import os
from string import Template

import pandas as pd
pd.set_option('display.max_colwidth', -1)

from IPython.core.display import HTML

from . import templ
from .mol_images import mol_img_tag

LIB_LOCATION = "local"  # other option: "net"
# for local to work, it requires the following files to be available
# in lib/css:
# jquery.dataTables.min.css
# bootstrap1.min.css
# bootstrap2.min.css
# in lib/:
# jquery.min.js
# popper.js
# bootstrap.min.js
# jquery.dataTables.min.js

# If the files are not found, the module defaults to loading the web versions,
# which is probably what you want anyway.

SHOW_WARN = True


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


def _mol_img_tag(mol):
    return pd.Series(mol_img_tag(mol))


def df_html(df, title="MolFrame", include_smiles=False,
            drop=[], keep=[], fn="tmp.html", **kwargs):
    df = df.copy()
    """Known kwargs: smiles_col, mol_col, id_col, b64_col, fp_col (str); index (bool)"""
    # ---- KW Args ----
    smiles_col = kwargs.get("smiles_col", "Smiles")
    mol_col = kwargs.get("mol_col", "Mol")
    id_col = kwargs.get("id_col", "Compound_Id")
    b64_col = kwargs.get("b64_col", "Mol_b64")
    fp_col = kwargs.get("fp_col", "FP_b64")
    index = kwargs.get("index", False)
    drop.extend([b64_col, fp_col])
    if not include_smiles:
        drop.append(smiles_col)
    if len(keep) > 0:
        keep.append(mol_col)
        df = df[keep]
    keys = list(df.keys())
    mol_col_pos = keys.index(mol_col)
    keys.pop(mol_col_pos)
    keys_sort = [mol_col]
    if id_col in keys:
        id_col_pos = keys.index(id_col)
        keys_sort.append(id_col)
        keys.pop(id_col_pos)
    keys_sort.extend(keys)
    df = df[keys_sort]
    if len(drop) > 0:
        # Find the keys to drop that are actually still in the df
        drop = list(set(drop).intersection(set(df.keys())))
        df = df.drop(drop, axis=1)
    tbl = df.to_html(formatters={mol_col: _mol_img_tag}, escape=False, index=index)
    tbl = tbl.replace("<td>0    <img", "<td><img")
    tbl = tbl.replace("dtype: object", "")
    return tbl


def view(df, title="MolFrame", include_smiles=False, drop=[], keep=[], fn="tmp.html", **kwargs):
    """Known kwargs: smiles_col, mol_col, id_col, b64_col (str); selectable (bool), index (bool),
        intro (text)"""
    global SHOW_WARN
    if "local" in LIB_LOCATION.lower() and os.access("lib/bootstrap.min.js", os.R_OK):
        pandas_tbl = templ.PANDAS_TABLE_LOCAL
    else:
        pandas_tbl = templ.PANDAS_TABLE_NET
        if SHOW_WARN:
            SHOW_WARN = False
            print("* using online libs for dataframe browsing...")
    intro = kwargs.get("intro", "")
    tbl = df_html(df, title, include_smiles, drop, keep, fn, **kwargs)
    tbl_list = tbl.split("\n")
    tbl_list = tbl_list[1:-1]
    tbl = "\n".join(tbl_list)
    templ_dict = {"title": title, "table": tbl, "intro": intro}

    if kwargs.get("selectable", False):
        templ_dict["selection_btn"] = templ.SELECTION_BTN
        templ_dict["selection_js"] = templ.SELECTION_JS
    else:
        templ_dict["selection_btn"] = ""
        templ_dict["selection_js"] = ""

    t = Template(pandas_tbl)
    html = t.substitute(templ_dict)
    write(html, fn)
    return HTML('<a href="{}">{}</a>'.format(fn, title))
