#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
##############
Resource Paths
##############

*Created on Sun Aug 6, 2017 18:30 by A. Pahl*

This is a template file. Adapt to your needs and rename to resource_paths.py
"""

# tab-delim., gzipped
smiles_path = "xxx/smiles_b64.tsv.gz"
smiles_cols = ['Compound_Id', "Smiles"]

# Batch based data, like purity, producer; tab-delim., gzipped
batch_data_path = "xxx/batch_data_b64.tsv.gz"
batch_data_cols = ["Compound_Id", 'Batch_Id', "Producer", "Pure_Flag", "Avail"]

# container based data, e.g. biological activities; tab-delim., gzipped
container_data_path = "xxx/container_data_b64.tsv.gz"
container_data_cols = []
