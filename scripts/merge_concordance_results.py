#!/usr/bin/env python3

"""
Merges the per-chromosome results of one threshold combination: the confusion matrices, the variant
missing-rate tables (vmiss) and the sample missing-rate tables (smiss) stored in several .pkl files. From each
file it takes the entries whose thresholds (DP, GQ, LAF, HAF) match the command-line arguments, and merges them
separately for each sample group (ALL, 15X, 30X).

How the tables are merged:
- vmiss: concatenated and sorted by VARIANT_ID, first by chromosome (chr1-chr22), then by position.
- smiss: grouped by SAMPLE_ID; SMISS_COUNT and TOTAL_VARIANT are summed and SMISS_FREQ recomputed.
- confusion: the counts of identical (CALL_GENOTYPE, TRUE_GENOTYPE) pairs are added.

Arguments:
--dp            DP threshold (int)
--gq            GQ threshold (int)
--laf           low allele-fraction bound, LAF (float)
--haf           high allele-fraction bound, HAF (float)
--input_chr_pkl the .pkl files, one per chromosome

Output:
- one .pkl file named merged._DP{dp}_GQ{gq}_LAF{laf}_HAF{haf}_.pkl, holding for each sample group (ALL, 15X,
  30X) the merged vmiss, smiss and confusion tables.

Example:
python3 merge_concordance_results.py --dp 10 --gq 20 --laf 0.1 --haf 0.9 --input_chr_pkl chr1.pkl chr2.pkl ...
"""

import argparse
import pickle
import pandas as pd
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(description="Merge confusion, vmiss and smiss data from multiple .pkl files")
parser.add_argument('--dp', type=int, required=True)
parser.add_argument('--gq', type=int, required=True)
parser.add_argument('--laf', type=float, required=True)
parser.add_argument('--haf', type=float, required=True)
parser.add_argument('--input_chr_pkl', nargs='+', required=True, help='List of input .pkl files (one per chromosome)')
args = parser.parse_args()

def merge_tag(tag, args, all_vmiss, all_smiss, all_confusion):
    if all_vmiss[tag]:
        merged_vmiss_df = pd.concat(all_vmiss[tag], ignore_index=True)
        # sort by chromosome, then position
        if not merged_vmiss_df.empty and "VARIANT_ID" in merged_vmiss_df.columns:
            merged_vmiss_df["__CHROM_INT"] = merged_vmiss_df["VARIANT_ID"].str.split(":").str[0].str.replace("chr", "").astype(int)
            merged_vmiss_df["__POS_INT"] = merged_vmiss_df["VARIANT_ID"].str.split(":").str[1].astype(int)
            merged_vmiss_df.sort_values(["__CHROM_INT", "__POS_INT"], inplace=True)
            merged_vmiss_df.drop(columns=["__CHROM_INT", "__POS_INT"], inplace=True)
    else:
        merged_vmiss_df = pd.DataFrame()

    if all_smiss[tag]:
        smiss_merged_df = pd.concat(all_smiss[tag], ignore_index=True)
        if not smiss_merged_df.empty:
            smiss_agg_df = (
                smiss_merged_df
                .groupby("SAMPLE_ID", as_index=False)
                .agg({
                    "SMISS_COUNT": "sum",
                    "TOTAL_VARIANT": "sum"
                })
            )
            smiss_agg_df["SMISS_FREQ"] = smiss_agg_df["SMISS_COUNT"] / smiss_agg_df["TOTAL_VARIANT"]
        else:
            smiss_agg_df = pd.DataFrame()
    else:
        smiss_agg_df = pd.DataFrame()

    confusion_df_combined = pd.DataFrame([
        {"CALL_GENOTYPE": k[0], "TRUE_GENOTYPE": k[1], "COUNT": v}
        for k, v in all_confusion[tag].items()
    ])

    output_key = frozenset([
        f'DP{args.dp}',
        f'GQ{args.gq}',
        f'LAF{args.laf}',
        f'HAF{args.haf}',
        tag
    ])
    return output_key, (merged_vmiss_df, smiss_agg_df, confusion_df_combined)



tags = ["ALL", "15X", "30X"]
param_keys = {
    tag: frozenset([
        f'DP{args.dp}',
        f'GQ{args.gq}',
        f'LAF{args.laf}',
        f'HAF{args.haf}',
        tag
    ]) for tag in tags
}

# === Containers ===
all_vmiss = {tag: [] for tag in tags}
all_smiss = {tag: [] for tag in tags}
all_confusion = {tag: defaultdict(int) for tag in tags}

# === Read each .pkl file ===
for pkl_file in args.input_chr_pkl:
    with open(pkl_file, "rb") as f:
        data = pickle.load(f)

    for key, (vmiss_df, smiss_df, confusion_df) in data.items():
        for tag, param_key in param_keys.items():
            if key == param_key:
                if vmiss_df is not None:
                    all_vmiss[tag].append(vmiss_df)
                if smiss_df is not None:
                    all_smiss[tag].append(smiss_df)
                if confusion_df is not None:
                    for _, row in confusion_df.iterrows():
                        tup = (row["CALL_GENOTYPE"], row["TRUE_GENOTYPE"])
                        all_confusion[tag][tup] += row["COUNT"]

with ProcessPoolExecutor() as executor:
    futures = [executor.submit(merge_tag, tag, args, all_vmiss, all_smiss, all_confusion) for tag in tags]
    merged_results_dict = dict(f.result() for f in futures)

output_path = f"merged._DP{args.dp}_GQ{args.gq}_LAF{args.laf}_HAF{args.haf}_.pkl"
with open(output_path, "wb") as out_f:
    pickle.dump(merged_results_dict, out_f)
print(f"[✓] Merged results for ALL/15X/30X saved → {output_path}")