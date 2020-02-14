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

import os.path as op
import time

# from string import Template

import pandas as pd

pd.set_option("display.max_colwidth", None)


from mol_frame import templ
from mol_frame import nb_tools as nbt
from mol_frame.mol_images import b64_mol, mol_img_tag, mol_img_file

_ = mol_img_file

IPYTHON = nbt.is_interactive_ipython()
if IPYTHON:
    from IPython.core.display import HTML

BGCOLOR = "#94CAEF"
IMG_GRID_SIZE = 235

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


def _apply_link(input, link, ln_title="Link"):
    """input[0]: mol_img_tag
    input[1]: link_col value"""
    link_str = link.format(input[1])
    result = '<a target="_blank" href="{}" title="{}">{}</a>'.format(
        link_str, ln_title, input[0]
    )
    return result


def drop_cols(df, drop):
    """Drop only exiosting cols from df.
    Raises no error, when cols do not exist."""
    if len(drop) > 0:
        # Find the keys to drop that are actually still in the df
        drop = list(set(drop).intersection(set(df.keys())))
        if len(drop) > 0:
            df = df.drop(drop, axis=1)
    return df


def show_mols(mols_or_smiles, cols=4):
    """A small utility to quickly view a list of mols in a grid."""
    html_list = ['<table align="center">']
    idx = 0
    row_line = []
    if not isinstance(mols_or_smiles, list):
        mols_or_smiles = [mols_or_smiles]
    for mol in mols_or_smiles:
        idx += 1
        img = mol_img_tag(mol)
        cell = "<td>{}<td>".format(img)
        row_line.append(cell)
        if idx == cols:
            row = "<tr>" + "".join(row_line) + "</td>"
            html_list.append(row)
            idx = 0
            row_line = []
    if idx != 0:
        row = "<tr>" + "".join(row_line) + "</td>"
        html_list.append(row)
    html_list.append("</table>")
    table = "\n".join(html_list)
    return HTML(table)


def html_grid(
    df,
    title="MolFrame",
    drop=[],
    keep=[],
    smiles_col="Smiles",
    mol_col="Mol",
    id_col=None,
    b64_col=None,
    fp_col=None,
    size=IMG_GRID_SIZE,
    **kwargs,
):
    """Creates a HTML grid out of the MolFrame.data input.

    Parameters:
        df (MolFrame): list of RDKit molecules
        highlight (dict): dict of properties (a.t.m only one) and values to highlight cells,
            e.g. {"activity": "< 50"}
        order (list): a list of substrings to match with the field names for ordering in the table header
        img_dir (str): if None, the molecule images are embedded in the HTML doc.
            Otherwise the images will be stored in img_dir and linked in the doc.
        interactive (bool)
        link_templ, link_col (str) (then interactive is false)
    Returns:
        HTML table as TEXT with molecules in grid-like layout to embed in IPython or a web page."""
    interact = kwargs.get("interactive", False)
    link_templ = kwargs.get("link_templ", None)
    link_col = kwargs.get("link_col", None)
    if link_col is not None:
        interact = (
            False  # interact is not avail. when clicking the image should open a link
        )
    highlight = kwargs.get("highlight", None)
    mols_per_row = kwargs.get("mols_per_row", 4)
    hlsss = kwargs.get(
        "hlsss", None
    )  # colname with Smiles (,-separated) for Atom highlighting
    truncate = kwargs.get("truncate", 12)
    img_dir = None
    if len(keep) > 0:
        keep.append(mol_col)
        if id_col is not None and id_col not in keep:
            keep.append(id_col)
        df = df[keep]
    drop.extend([smiles_col, b64_col, fp_col])
    df = drop_cols(df, drop)
    props = []
    for x in list(df.keys()):
        if x != mol_col and x != id_col and x != hlsss:
            props.append(x)
    time_stamp = time.strftime("%y%m%d%H%M%S")
    # td_opt = {"align": "center"}
    td_opt = {"style": "text-align: center;"}
    header_opt = {"bgcolor": BGCOLOR}
    table_list = []
    guessed_id = id_col
    if interact and guessed_id is not None:
        table_list.append(templ.TBL_JAVASCRIPT.format(ts=time_stamp, bgcolor=BGCOLOR))

    if len(props) > 0:
        td_opt["colspan"] = "2"
        prop_row_cells = {k: [] for k, _ in enumerate(props)}

    rows = []
    id_cells = []
    mol_cells = []
    for idx, (_, rec) in enumerate(df.iterrows(), 1):
        mol = rec[mol_col]
        if guessed_id:
            id_prop_val = str(rec[guessed_id])
            img_id = id_prop_val
            cell_opt = {"id": "{}_{}".format(id_prop_val, time_stamp)}
            cell_opt.update(td_opt)
            cell_opt.update(header_opt)
            id_cells.extend(templ.td(id_prop_val, cell_opt))
        else:
            img_id = idx

        if not mol:
            cell = ["no structure"]

        else:
            if hlsss is not None:
                hlsss_smi = rec[hlsss]
            else:
                hlsss_smi = None
            if img_dir is None:  # embed the images in the doc
                img_src = b64_mol(mol, size * 2, hlsss=hlsss_smi)

            if interact and guessed_id is not None:
                img_opt = {
                    "title": "Click to select / unselect",
                    "onclick": "toggleCpd('{}')".format(id_prop_val),
                }
            elif link_col is not None:

                img_opt = {"title": "Click to open link"}
                #          "onclick": "location.href='{}';".format(link)}
                # '<a target="_blank" href="{}" title="{}">{}</a>'

            else:
                img_opt = {"title": str(img_id)}

            img_opt["style"] = f"max-width: {size}px; max-height: {size}px;"
            # img_opt["height"] = "{}px".format(size)
            cell = templ.img(img_src, img_opt)
            if link_col is not None:
                link = link_templ.format(rec[link_col])
                a_opt = {"href": link}
                cell = templ.a(cell, a_opt)

        # td_opt = {"align": "center"}
        td_opt = {"style": "text-align: center;", "bgcolor": "#FFFFFF"}
        if len(props) > 0:
            td_opt["colspan"] = "2"

        if highlight is not None:
            eval_str = None
            prop = highlight.keys()[0]  # only one highlight key supported a.t.m.
            prop_val = str(rec[prop])
            eval_str = " ".join([prop_val, highlight[prop]])
            if eval_str and eval(eval_str):
                td_opt["bgcolor"] = "#99ff99"

        mol_cells.extend(templ.td(cell, td_opt))

        if len(props) > 0:
            for prop_no, prop in enumerate(props):
                prop_opt = {"style": "text-align: left;"}
                val_opt = {"style": "text-align: left;"}
                prop_cells = []
                prop_val = ""
                if prop in rec:
                    prop_val = str(rec[prop])
                    if (
                        prop == "Pure_Flag"
                        and prop_val != ""
                        and prop_val != "n.d."
                        and "Purity" in rec
                        and "LCMS_Date" in rec
                    ):
                        val_opt["title"] = "{}% ({})".format(
                            rec["Purity"], rec["LCMS_Date"]
                        )
                prop_cells.extend(templ.td(prop[:25], prop_opt))
                prop_cells.extend(templ.td(prop_val[:truncate], val_opt))
                prop_row_cells[prop_no].extend(prop_cells)

        if idx % mols_per_row == 0 or idx == len(df):
            if guessed_id:
                rows.extend(templ.tr(id_cells))
            rows.extend(templ.tr(mol_cells))

            if len(props) > 0:
                colspan_factor = 2
                for prop_no in sorted(prop_row_cells):
                    rows.extend(templ.tr(prop_row_cells[prop_no]))
                prop_row_cells = {k: [] for k, _ in enumerate(props)}
            else:
                colspan_factor = 1
            empty_row_options = {"colspan": mols_per_row * colspan_factor}
            empty_row_options["style"] = "border: none;"
            empty_row = templ.tr(templ.td("&nbsp;", options=empty_row_options))
            rows.extend(empty_row)
            id_cells = []
            mol_cells = []

    table_list.extend(templ.table(rows))

    if interact and guessed_id is not None:
        table_list.append(templ.ID_LIST.format(ts=time_stamp))

    # print(table_list)
    return "".join(table_list)


def write_grid(
    df,
    title="MolGrid",
    fn="molgrid.html",
    drop=[],
    keep=[],
    smiles_col="Smiles",
    mol_col="Mol",
    id_col=None,
    b64_col=None,
    fp_col=None,
    **kwargs,
):
    """Write the html_grid to file and return the link."""
    tbl = html_grid(
        df,
        title=title,
        drop=drop,
        keep=keep,
        smiles_col=smiles_col,
        mol_col=mol_col,
        id_col=id_col,
        b64_col=b64_col,
        fp_col=fp_col,
        size=IMG_GRID_SIZE,
        **kwargs,
    )
    header = kwargs.get("header", None)
    summary = kwargs.get("summary", None)
    page = templ.page(tbl, title=title, header=header, summary=summary)
    write(page, fn)
    if IPYTHON:
        return HTML('<a href="{}">{}</a>'.format(fn, title))


def show_grid(
    df,
    title="MolFrame",
    drop=[],
    keep=[],
    smiles_col="Smiles",
    mol_col="Mol",
    id_col=None,
    b64_col=None,
    fp_col=None,
    **kwargs,
):
    """Show the html_grid in the Notebook."""
    html = html_grid(
        df,
        title=title,
        drop=drop,
        keep=keep,
        smiles_col=smiles_col,
        mol_col=mol_col,
        id_col=id_col,
        b64_col=b64_col,
        fp_col=fp_col,
        size=IMG_GRID_SIZE,
        **kwargs,
    )
    if IPYTHON:
        return HTML(html)
    else:
        return html


def rm_table_tag(tbl):
    tbl_list = tbl.split("\n")
    tbl_list = tbl_list[1:-1]
    result = "\n".join(tbl_list)
    return result


def jsme(name="mol"):
    """Displays a JSME molecule editor widget in the notebook
    and stores the resulting mol in the variable that <name> assigns.
    Requires the following line to be imported in the Notebook:
    ``from rdkit.Chem import AllChem as Chem``"""
    if not IPYTHON:
        print("* interactive Jupyter notebook session required.")
        return
    if op.isfile("lib/jsme/jsme.nocache.js"):
        JSME_LOCATION = "lib"
    else:
        print("* no local installation of JSME found, using web version.")
        JSME_LOCATION = "http://peter-ertl.com/jsme/JSME_2017-02-26"

    time_stamp = time.strftime("%y%m%d%H%M%S")

    return HTML(
        templ.JSME_FORM.format(jsme_loc=JSME_LOCATION, ts=time_stamp, var_name=name)
    )
