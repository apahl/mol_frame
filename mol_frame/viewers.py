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
import base64
from io import BytesIO as IO

import pandas as pd
pd.set_option('display.max_colwidth', -1)
from PIL import Image, ImageChops

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
Draw.DrawingOptions.atomLabelFontSize = 18

from IPython.core.display import HTML

from . import templ

LIB_LOCATION = "local"
SHOW_WARN = True
print("* viewers reloaded")


def write(text, fn):
    with open(fn, "w") as f:
        f.write(text)


def make_transparent(img):
    img = img.convert("RGBA")
    pixdata = img.load()
    width, height = img.size
    for y in range(height):
        for x in range(width):
            if pixdata[x, y] == (255, 255, 255, 255):
                pixdata[x, y] = (255, 255, 255, 0)
    return img


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    return None  # no contents


def b64_mol(mol, size=300):
    img_file = IO()
    try:
        img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    except UnicodeEncodeError:
        print(Chem.MolToSmiles(mol))
        mol = Chem.MolFromSmiles("*")
        img = autocrop(Draw.MolToImage(mol, size=(size, size)))
    img = make_transparent(img)
    img.save(img_file, format='PNG')
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def mol_img_tag(mol, options=None):
    tag = """<img {} src="data:image/png;base64,{}" alt="Mol"/>"""
    if options is None:
        options = ""
    img_tag = tag.format(options, b64_mol(mol))
    return img_tag


def _mol_img_tag(mol):
    return pd.Series(mol_img_tag(mol))


def df_html(df, title="MolFrame", include_smiles=False,
            drop=[], keep=[], fn="tmp.html", **kwargs):
    df = df.copy()
    """Known kwargs: smiles_col, mol_col, id_col"""
    # ---- KW Args ----
    smiles_col = kwargs.get("smiles_col", "Smiles")
    mol_col = kwargs.get("mol_col", "Mol")
    id_col = kwargs.get("id_col", "Compound_Id")
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
    tbl = df.to_html(formatters={mol_col: _mol_img_tag}, escape=False)
    tbl = tbl.replace("<td>0    <img", "<td><img")
    tbl = tbl.replace("dtype: object", "")
    return tbl


def view(df, title="MolFrame", include_smiles=False, drop=[], keep=[], fn="tmp.html", **kwargs):
    """Known kwargs: smiles_col, mol_col, id_col"""
    global SHOW_WARN
    if "local" in LIB_LOCATION.lower() and os.access("lib/bootstrap.min.js", os.R_OK):
        pandas_tbl = templ.PANDAS_TABLE_LOCAL
    else:
        pandas_tbl = templ.PANDAS_TABLE_NET
        if SHOW_WARN:
            SHOW_WARN = False
            print("* using online libs for dataframe browsing...")
    tbl = df_html(df, title, include_smiles, drop, keep, fn, **kwargs)
    tbl_list = tbl.split("\n")
    tbl_list = tbl_list[1:-1]
    tbl = "\n".join(tbl_list)
    t = Template(pandas_tbl)
    html = t.substitute(title=title, table=tbl)
    write(html, fn)
    return HTML('<a href="{}">{}</a>'.format(fn, title))
