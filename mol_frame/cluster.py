#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
#######
Cluster
#######

*Created on Thu Mar 28, 2019 by A. Pahl*

Tools for clustering structures."""

import pandas as pd

from mol_frame import mol_frame as mf, templ as mfht, viewers as mfv

from typing import List


def report(
    df: pd.DataFrame,
    id_col: str = "Compound_Id",
    columns: List[str] = ["Compound_Id", "Smiles"],
    title: str = "Cluster Report",
    intro: str = "Large clusters first, similar clusters together.",
):
    """Write a HTML report. `Cluster_No` and `IsRepr` have to be present in the DataFrame.
    In the current setting, the largest clusters are at the top of the report,
    with similar clusters (determind by the chemical similarities of the representative structures)
    are grouped together.
    Writes the report to disk as `Clusters.html`.
    Used in `projects/paint3_anal/190328_cpd_clustering.ipynb`.

    Arguments:
        df: The input DataFrame containing the structures as Smiles.
        id_col: The name of the column to use for identity. Default is `Compound_Id`.
        columns: List of columns to include.
        title: The report title.
        intro: Some text used for introduction of the report.
    """

    def add_cluster(cl_no, sim_to=None):
        if sim_to is None:
            sim_to = ""
            html.append("<hr>")
        else:
            sim_to = f"(similar to {sim_to})"
        mf_cl = mf.MolFrame(df.query("Cluster_No == @cl_no")[columns])
        mf_cl = mf_cl.add_mols()
        html.append(
            f"<br><h2>Cluster {cl_no} ({len(mf_cl.data)} Members)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;{sim_to}</h2><br>"
        )
        grid = mfv.html_grid(mf_cl.data, id_col="Compound_Id")
        html.append(grid)

    if id_col not in columns:
        columns = [id_col] + columns
    if "Smiles" not in columns:
        columns.append("Smiles")
    df_repr = df.query("IsRepr == 'Yes'").reset_index().drop("index", axis=1)
    chem_sim = {}
    for idx, rec0 in df_repr.iterrows():
        for _, rec1 in df_repr.iloc[idx + 1 :].iterrows():
            cl0 = rec0["Cluster_No"]
            cl1 = rec1["Cluster_No"]
            sim = mf.chem_sim(rec0["Smiles"], rec1["Smiles"])
            chem_sim[(cl0, cl1)] = sim
            chem_sim[(cl1, cl0)] = sim

    cl_sizes = (
        df[["Cluster_No", "Compound_Id"]]
        .groupby(by="Cluster_No")
        .count()
        .reset_index()
        .rename(columns={"Compound_Id": "Size"})
    )
    cl_sizes = cl_sizes.sort_values("Size", ascending=False)
    cl_order = {x: True for x in cl_sizes["Cluster_No"].values}

    html = [f"<h1>{title}</h1><br>{intro}<br><br>"]
    while len(cl_order) > 0:
        cl_no = list(cl_order.keys())[0]
        add_cluster(cl_no)
        cl_order.pop(cl_no)
        to_remove = []
        for sim_cl in cl_order:
            if chem_sim[(cl_no, sim_cl)] > 0.45:
                add_cluster(sim_cl, cl_no)
                to_remove.append(sim_cl)
        for x in to_remove:
            cl_order.pop(x)

    mfht.write(mfht.page("\n".join(html)), "Clusters.html")
