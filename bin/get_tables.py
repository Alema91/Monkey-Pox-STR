#!/usr/bin/env python

# Imports
import sys

import pandas as pd
from sklearn import preprocessing


# Arguments
raw_table = sys.argv[1]

# load table, change sequence label for a number
raw_df = pd.read_excel(
    raw_table, 
    index_col = "Sample_name")

le = preprocessing.LabelEncoder()
raw_df["Sequence"] = le.fit_transform(raw_df["Sequence"])

# 1: only highest frequency alleles

# List of all STRs - Samples
STRS = list(set(raw_df["STR_mark"]))
SAMPLES = list(set(raw_df.index))

# dicts for allele freq
alfreq_dict = dict()
seq_dict = dict()

for STR in STRS:
    alfreq_dict[STR] = dict()
    seq_dict[STR] = dict()
    for sample in SAMPLES:
        alfreq_dict[STR][sample] = 0
        seq_dict[STR][sample] = "None"

# find most frequent alleles, crate df with them
for index, row in raw_df.iterrows():
    STR = row["STR_mark"]
    al_freq = row["AlleleFrequency"]
    sample = index
    if alfreq_dict[STR][sample] < al_freq:
        seq = row["Sequence"]
        alfreq_dict[STR][sample] = al_freq
        seq_dict[STR][sample] = seq

df_xsample_ySTR = pd.DataFrame.from_dict(seq_dict)
df_xsample_ySTR.to_csv("01-Most_abundant_alelles.tsv", sep="\t")
