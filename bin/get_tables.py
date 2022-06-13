#!/usr/bin/env python

# Imports
import os
import sys
import json
import pandas as pd
from composition_stats import closure
from composition_stats import clr


# Required functions
def get_STRS_samplelist(df):
    # Get STR and sample list
    STRS = list(set(df["STR_mark"]))
    sample_list = list(set(df.index))
    
    return STRS, sample_list

def encode_sequences(df):
    traduction_dict = dict()
    count_dict = dict()
    
    for row in range(df.shape[0]):

        STR = df.iloc[row]["STR_mark"]
        seq = df.iloc[row]["Sequence"]

        if STR not in traduction_dict.keys():
            traduction_dict[STR] = {}
            count_dict[STR] = 0 

        if seq not in traduction_dict[STR].keys():
            new_seq = f"{STR}_allele_{count_dict[STR]}"
            traduction_dict[STR][seq] = new_seq
            # Sequence = col 5
            df.iloc[row,5] = new_seq

            count_dict[STR]+=1

        else:

            df.iloc[row,5] = traduction_dict[STR][seq]
        
    return df, traduction_dict
    

def most_frequent_alleles(df):

    STRS, sample_list = get_STRS_samplelist(df)
    
    # alfreq_dict, relfreqs of alleles, starts at 0
    # seq_dict, seq of alleles, starts at None
    alfreq_dict = dict()
    seq_dict = dict()

    # get the most frequent allele and its relfreq
    for STR in STRS:
        alfreq_dict[STR] = dict()
        seq_dict[STR] = dict()
        for sample in sample_list:
            alfreq_dict[STR][sample] = 0
            seq_dict[STR][sample] = "None"
    
    for index, row in df.iterrows():
        STR = row["STR_mark"]
        al_freq = row["AlleleFrequency"]
        sample = index

        if alfreq_dict[STR][sample] < al_freq:
            seq = row["Sequence"]

            alfreq_dict[STR][sample] = al_freq
            seq_dict[STR][sample] = seq
    
    # seq of the most freq allele per str and sample
    df_mostfreqalleles_seqs = pd.DataFrame.from_dict(seq_dict)
    
    # freq of the most freq allele per str and sample
    df_mostfreqalleles_freqs = pd.DataFrame.from_dict(alfreq_dict)

            
    return df_mostfreqalleles_seqs, df_mostfreqalleles_freqs
    
def all_alleles(df, column):
    STRS, sample_list = get_STRS_samplelist(df)
    
    # rename alleles
    raw_seq_renamed_df = raw_df.copy()
    
    # list with all different alleles
    all_alleles = list(set(raw_seq_renamed_df["Sequence"]))
    
    # dict for freq of every allele
    alleles_dict_freq = { item: { allele:0 for allele in all_alleles } for item in sample_list }
    alleles_dict_binary = { item: { allele:0 for allele in all_alleles } for item in sample_list }

    # dict for presence of every allele
    
    for index, row in raw_seq_renamed_df.iterrows():
        sample = index
        allele = row["Sequence"]
        freq = row[column]

        alleles_dict_freq[sample][allele] = float(freq)
        alleles_dict_binary[sample][allele] = 1
    
        # df for presence
    df_alleles_binary = pd.DataFrame.from_dict(alleles_dict_binary)
    
    # df for freq
    df_alleles_freq = pd.DataFrame.from_dict(alleles_dict_freq)
        
    return df_alleles_binary, df_alleles_freq

def clr_df(df, pseudosum=0.000001):
    rows = list(df.columns)
    # alleles will be the cols
    cols = list(df.index)
    clr_df = pd.DataFrame(clr(closure(df_alleles_supporting_reads.transpose()+pseudosum)),columns=cols, index=rows)

    return clr_df

# Arguments
raw_table = sys.argv[1]

# Import data
raw_df = pd.read_excel(raw_table, index_col = "Sample_name")

# Encode data
renamed_df, traduction_dict = encode_sequences(raw_df)

# save enconding dict as json
with open("Encoding_dict.json", "w") as outfile:
    json.dump(traduction_dict, outfile)

# tables for most frequent alleles
df_mostfreqalleles_seqs, df_mostfreqalleles_freqs = most_frequent_alleles(renamed_df)
df_mostfreqalleles_seqs.to_csv("01-most_frequent_alleles_seqs.tsv", sep = "\t")
df_mostfreqalleles_freqs.to_csv("02-most_frequent_alleles_freqs.tsv")

# tables for all alleles
df_alleles_binary, df_alleles_freq = all_alleles(renamed_df, "AlleleFrequency")
df_alleles_binary.to_csv("03-all_alleles_presence.tsv", sep="\t")
df_alleles_freq.to_csv("04-all_alleles_relfreq.tsv", sep="\t")

# table for all alleles by supporting reads
df_alleles_supporting_reads = all_alleles(renamed_df, "Supporting_reads")[1]
df_alleles_supporting_reads.to_csv("05-all_alleles_supporting_reads.tsv", sep="\t")

# apply a clr conversion to supporting reads
df_clr = clr_df(df_alleles_supporting_reads)
df_clr.to_csv("06-all_alleles_supporting_reads_CLR.tsv", sep="\t")

