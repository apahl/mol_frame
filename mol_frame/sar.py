#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###
SAR
###

*Created on Thu Mar 28, 2019 by A. Pahl*

Tools for SAR analysis."""

import base64, pickle, time
from io import BytesIO as IO
import os.path as op
from collections import Counter
from copy import deepcopy

# import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit import DataStructs

try:
    Draw.DrawingOptions.atomLabelFontFace = "DejaVu Sans"
    Draw.DrawingOptions.atomLabelFontSize = 18
except KeyError:  # Font "DejaVu Sans" is not available
    pass

from mol_frame import mol_frame as mf, mol_images as mi

# from typing import List,


class SAR:
    """Container class for SAR analysis.
    Operates on a copy of the original MolFrame.
    All methods of the SAR class return copies,
    except for the `train` method."""

    def __init__(self, molf: mf.MolFrame = None):
        """
        Parameters:
            molf: MolFrame instance."""

        if molf is not None:
            self.molf = molf.copy()
        self.model = None

    def __str__(self):
        shape = self.molf.data.shape
        keys = list(self.molf.data.keys())
        return f"MolFrame  Rows: {shape[0]:6d}  Columns: {shape[1]:2d}   {keys}"

    def __repr__(self):
        return self.__str__()

    def __getitem__(self, key):
        res = self.molf[key]
        if isinstance(res, mf.MolFrame):
            result = self.new()
            result.molf = res
        else:
            result = res
        return result

    def __getattr__(self, name):
        """Try to call undefined methods on the underlying pandas DataFrame."""
        if hasattr(self.molf, name):

            def method(*args, **kwargs):
                res = getattr(self.molf, name)(*args, **kwargs)
                if isinstance(res, mf.MolFrame):
                    result = self.new()
                    result.molf = res
                else:
                    result = res
                return result

            return method
        else:
            raise AttributeError

    def new(self):
        result = SAR()
        if self.model is None:
            result.model = None
        else:
            result.model = deepcopy(self.model)
        return result

    def write(self, **kwargs):
        bn = kwargs.get("name", self.config["NAME"])
        self.molf.write_csv(f"{bn}.tsv")

    def copy(self):
        result = SAR(self.molf)
        result.model = self.model
        return result

    def to_csv(self, fn, sep="\t"):
        self.molf.write_csv(fn, sep="\t")

    def analyze(self, act_class="AC_Real", pred_class="AC_Pred"):
        """Prints the ratio of succcessful predictions for the molecules which have `act_class` and `pred_class` properties."""
        mol_ctr = Counter()
        hit_ctr = Counter()
        for _, rec in self.molf.data.iterrows():
            if act_class in rec and pred_class in rec:
                mol_ctr[rec[act_class]] += 1
                if rec[act_class] != rec[pred_class]:
                    continue
                hit_ctr[rec[act_class]] += 1
        if len(mol_ctr) > 0:
            sum_mol_ctr = sum(mol_ctr.values())
            sum_hit_ctr = sum(hit_ctr.values())
            print(
                "Number of correctly predicted molecules: {} / {}    ({:.2f}%)".format(
                    sum_hit_ctr, sum_mol_ctr, 100 * sum_hit_ctr / sum_mol_ctr
                )
            )
            print("\nCorrectly predicted molecules per Activity Class:")
            for c in sorted(hit_ctr):
                print("  {}:  {:.2f}".format(c, 100 * hit_ctr[c] / mol_ctr[c]))
        else:
            print(
                "No molecules found with both {} and {}.".format(act_class, pred_class)
            )
        return hit_ctr, mol_ctr

    def save_model(self, fn="sar"):
        if self.model is None:
            print("No model available.")
            return
        save_model(self.model, fn)

    def load_model(self, fn="sar", force=False):
        if self.model is not None and not force:
            print("There is already a model available. Use `force=True` to override.")
            return
        if not fn.endswith(".model"):
            fn = fn + ".model"
        with open(fn, "rb") as f:
            self.model = pickle.load(f)
        print(
            "  > model loaded (last modified: {}).".format(
                time.strftime("%Y-%m-%d %H:%M", time.localtime(op.getmtime(fn)))
            )
        )

    def train(self, act_class="AC_Real", n_est=500, show_progress=True):
        self.model = train(
            self.molf, act_class=act_class, n_est=n_est, show_progress=show_progress
        )

    def predict(self):
        if self.model is None:
            raise LookupError("No suitable model found. Please run `train` first.")
        result = self.copy()
        result.molf = predict(self.molf, self.model)
        return result

    def add_sim_maps(self):
        """Adds the similarity maps as images to the MolFrame.
        Returns a copy."""

        result = self.copy()
        result.molf = add_sim_maps(self.molf, self.model)
        return result


def read_csv(name: str) -> SAR:
    bn = name
    molf = mf.read_csv(bn)
    result = SAR(molf)
    return result


def train(molf: mf.MolFrame, act_class="AC_Real", n_est=500, show_progress=True):
    """Returns the trained model."""
    fps = []
    act_classes = []
    molf.find_mol_col()
    if show_progress:
        print("  [TRAIN] calculating fingerprints")
    for _, rec in molf.data.iterrows():
        fps.append(
            Chem.GetMorganFingerprintAsBitVect(molf.mol_method(rec[molf.use_col]), 2)
        )
        act_classes.append(rec[act_class])
    np_fps = []
    if show_progress:
        print("  [TRAIN] calculating Numpy arrays")
    for fp in fps:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, arr)
        np_fps.append(arr)

    # get a random forest classifiert with 100 trees
    if show_progress:
        print("  [TRAIN] training RandomForestClassifier")
    rf = RandomForestClassifier(n_estimators=n_est, random_state=1123)
    rf.fit(np_fps, act_classes)
    if show_progress:
        print("  [TRAIN] done.")
    return rf


def predict(molf: mf.MolFrame, model):
    def _predict_mol(mol):
        """Returns the predicted class and the probabilities for a molecule.

        Parameters:
            model: Output from `train()`."""
        fp = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(Chem.GetMorganFingerprintAsBitVect(mol, 2), fp)
        fp = fp.reshape(1, -1)  # this removes the deprecation warning
        predict_class = model.predict(fp)
        predict_prob = model.predict_proba(fp)
        return predict_class[0], max(predict_prob[0])

    def _predict(s):
        mol = molf.mol_method(s[molf.use_col])
        result = _predict_mol(mol)
        return result  # returns tuple

    molf.find_mol_col()
    result = molf.copy()
    result.data[["AC_Pred", "Prob"]] = result.data.apply(
        _predict, axis=1, result_type="expand"
    )
    result.data["AC_Pred"] = result.data["AC_Pred"].astype(int)
    return result


def save_model(model, fn="sar"):
    if not fn.endswith(".model"):
        fn = fn + ".model"
    with open(fn, "wb") as f:
        pickle.dump(model, f)


def read_sdf(fn, model_name=None):
    sarf = SAR(mf.read_sdf(fn))
    if model_name is None:
        print("  * No model was loaded. Please provide a name to load.")
    else:
        try:
            sarf.load_model(model_name)
        except FileNotFoundError:
            print(
                "  * Model {} could not be found. No model was loaded".format(
                    model_name
                )
            )
    return sarf


def b64_fig(fig, dpi=72):
    img_file = IO()
    # print(fig.savefig.__doc__)
    # print([x for x in dir(fig) if x.startswith("set_")])
    # print(sorted(dir(fig)))
    # print(fig.get_edgecolor(), fig.get_facecolor())
    fig.savefig(img_file, dpi=dpi, format="PNG", bbox_inches="tight")
    img = mi.Image.open(img_file)
    img = mi.autocrop(img)
    img_file.close()
    img_file = IO()
    img.save(img_file, format="PNG")
    b64 = base64.b64encode(img_file.getvalue())
    b64 = b64.decode()
    img_file.close()
    return b64


def _get_proba(fp, predictionFunction):
    result = predictionFunction([fp])
    return result[0][1]


def add_sim_maps(molf: mf.MolFrame, model):
    """Adds the similarity maps as images to the MolFrame.
    Returns a copy."""

    def _map(x):
        mol = molf.mol_method(x)
        fig, _ = SimilarityMaps.GetSimilarityMapForModel(
            mol,
            SimilarityMaps.GetMorganFingerprint,
            lambda y: _get_proba(y, model.predict_proba),
            linewidths=0,
        )
        b64 = b64_fig(fig, dpi=72)
        img_src = '<img src="data:image/png;base64,{}" alt="Map" />'.format(b64)
        return img_src

    molf.find_mol_col()
    result = molf.copy()
    result.data["Map"] = result.data[molf.use_col].apply(_map)
    return result
