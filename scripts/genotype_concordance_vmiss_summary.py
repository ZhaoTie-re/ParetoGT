# %%
"""
Script: genotype concordance and missing-rate summary

Purpose:
Unpacks the evaluation results of one threshold combination (a .pkl from
evaluate_genotype_concordance_and_vmiss.py, or from merge_concordance_results.py) into three tables:
  - the variant missing-rate table (vmiss_df)
  - the sample missing-rate table (smiss_df)
  - the genotype confusion matrix (confusion_df)

and computes, for each sample group (ALL, 15X, 30X):
  - the number of genotypes
  - concordance, miss rate, false positive rate, false negative rate
  - the counts and rates of each kind of discordant call
  - the mean, standard deviation and median of VMISS and SMISS

Arguments:
  --chr: chromosome name (e.g. chr22), or 'merged' for the genome-wide result
  --global_dp: DP threshold
  --global_gq: GQ threshold
  --global_laf: low allele-fraction bound (LAF)
  --global_haf: high allele-fraction bound (HAF)
  --raw_pkl: the input .pkl, structured as
    {
      frozenset([...]): (vmiss_df, smiss_df, confusion_df)
    }

Output:
  - a .summary.pkl file, structured as
    {
      frozenset([...]): (vmiss_sub, smiss_sub, confusion_sub, summary_df)
    }
  with one entry per sample group.

"""
import pandas as pd
import numpy as np
import pickle
import argparse

# --- Arguments ---
parser = argparse.ArgumentParser(description="Genotype Concordance and VMISS Summary Script") 
parser.add_argument('--chr', type=str, required=True, help='Chromosome name (e.g., chr22)')
parser.add_argument('--global_dp', type=int, default=20, help='Global DP threshold')
parser.add_argument('--global_gq', type=int, default=30, help='Global GQ threshold')
parser.add_argument('--global_laf', type=float, default=0.25, help='Global LAF threshold')
parser.add_argument('--global_haf', type=float, default=0.75, help='Global HAF threshold')
parser.add_argument('--raw_pkl', type=str, required=True, help='Path to the raw .pkl file containing vmiss and confusion data')
args = parser.parse_args()  

# --- Settings ---
chr = args.chr
global_dp = args.global_dp
global_gq = args.global_gq
global_laf = args.global_laf
global_haf = args.global_haf

raw_pkl = args.raw_pkl
raw_dict = pd.read_pickle(raw_pkl)

# --- Unpacking ---
# Unpacks the result dictionary, adds the sample group and the thresholds as columns of each table, and
# concatenates the tables.
def unpack_concordance_vmiss_dict(raw_dict):
    all_vmiss, all_smiss, all_confusion = [], [], []

    for key, (vmiss_df, smiss_df, confusion_df) in raw_dict.items():
        key_strs = list(key)
        dp = int(next(s[2:] for s in key_strs if s.startswith('DP')))
        gq = int(next(s[2:] for s in key_strs if s.startswith('GQ')))
        laf = float(next(s[3:] for s in key_strs if s.startswith('LAF')))
        haf = float(next(s[3:] for s in key_strs if s.startswith('HAF')))
        platform = next(s for s in key_strs if s in ['ALL', '15X', '30X'])

        for df in (vmiss_df, smiss_df, confusion_df):
            df["DP"], df["GQ"], df["LAF"], df["HAF"], df["PLATFORM"] = dp, gq, laf, haf, platform

        all_vmiss.append(vmiss_df)
        all_smiss.append(smiss_df)
        all_confusion.append(confusion_df)

    return pd.concat(all_vmiss, ignore_index=True), pd.concat(all_smiss, ignore_index=True), pd.concat(all_confusion, ignore_index=True)

vmiss_df, smiss_df, confusion_df = unpack_concordance_vmiss_dict(raw_dict)

# --- Confusion matrix ---
# Turns confusion_df into a confusion matrix (rows: true genotype; columns: called genotype), with all four
# classes 0/1/2/NA present.
def confusion_to_matrix(confusion_df):
    ordered = ["0", "1", "2", "NA"]
    return (
        confusion_df.pivot_table(index="TRUE_GENOTYPE", columns="CALL_GENOTYPE", values="COUNT", aggfunc="sum", fill_value=0)
        .reindex(index=ordered, columns=ordered, fill_value=0)
    )

# --- Summary ---
# Computes the metrics from the confusion matrix and the missing-rate tables of one sample group.
import pandas as pd

def summarize_genotype_concordance(matrix, vmiss_df, smiss_df, global_dp, global_gq, global_laf, global_haf, platform):
    """
    Summarize genotype concordance from a confusion matrix and the missing-rate tables.

    Args:
        matrix: pd.DataFrame; rows TRUE_GENOTYPE, columns CALL_GENOTYPE, values COUNT
        global_dp, global_gq, global_laf, global_haf: the thresholds
        platform: str, the sample group, e.g. "15X", "30X", "ALL"
        vmiss_df: pd.DataFrame or None, with a "VMISS_FREQ" column
        smiss_df: pd.DataFrame or None, with a "SMISS_FREQ" column

    Returns:
        summary_df: a one-row DataFrame of the metrics
    """
    # index and columns as strings
    matrix = matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)

    # --- Totals ---
    true_values = ["0", "1", "2"]
    call_values = ["0", "1", "2"]

    total = matrix.loc[true_values, call_values].values.sum()
    total_with_empty = matrix.values.sum() - matrix.loc["NA"].sum() if "NA" in matrix.index else matrix.values.sum()

    # --- Concordance ---
    correct = sum(matrix.at[gt, gt] for gt in ["0", "1", "2"] if gt in matrix.index and gt in matrix.columns)
    concordance = correct / total if total > 0 else None
    concordance_with_empty = correct / total_with_empty if total_with_empty > 0 else None

    # --- Miss rate ---
    na_calls = matrix.loc[true_values, "NA"].sum() if "NA" in matrix.columns else 0
    miss_rate = na_calls / total_with_empty if total_with_empty > 0 else None

    # --- Discordant calls (false positives) ---
    het_to_homvar = matrix.at["1", "2"] if "1" in matrix.index and "2" in matrix.columns else 0
    homref_to_homvar = matrix.at["0", "2"] if "0" in matrix.index and "2" in matrix.columns else 0
    homref_to_het = matrix.at["0", "1"] if "0" in matrix.index and "1" in matrix.columns else 0

    het_to_homvar_rate = het_to_homvar / total if total > 0 else None
    homref_to_homvar_rate = homref_to_homvar / total if total > 0 else None
    homref_to_het_rate = homref_to_het / total if total > 0 else None
    
    # --- Discordant calls (false negatives) ---
    het_to_homref = matrix.at["1", "0"] if "1" in matrix.index and "0" in matrix.columns else 0
    homvar_to_homref = matrix.at["2", "0"] if "2" in matrix.index and "0" in matrix.columns else 0
    homvar_to_het = matrix.at["2", "1"] if "2" in matrix.index and "1" in matrix.columns else 0
    
    het_to_homref_rate = het_to_homref / total if total > 0 else None
    homvar_to_homref_rate = homvar_to_homref / total if total > 0 else None
    homvar_to_het_rate = homvar_to_het / total if total > 0 else None
    
    # --- Concordant calls ---
    homref_to_homref = matrix.at["0", "0"] if "0" in matrix.index and "0" in matrix.columns else 0
    het_to_het = matrix.at["1", "1"] if "1" in matrix.index and "1" in matrix.columns else 0
    homvar_to_homvar = matrix.at["2", "2"] if "2" in matrix.index and "2" in matrix.columns else 0
    
    homref_to_homref_rate = homref_to_homref / total if total > 0 else None
    het_to_het_rate = het_to_het / total if total > 0 else None
    homvar_to_homvar_rate = homvar_to_homvar / total if total > 0 else None

    # --- False positives and false negatives ---
    # False positive: the call carries more alternate alleles than the truth (0>1, 0>2, 1>2)
    false_positive = (
        matrix.at["0", "1"] if ("0" in matrix.index and "1" in matrix.columns) else 0
    ) + (
        matrix.at["0", "2"] if ("0" in matrix.index and "2" in matrix.columns) else 0
    ) + (
        matrix.at["1", "2"] if ("1" in matrix.index and "2" in matrix.columns) else 0
    )

    # False negative: the call carries fewer alternate alleles than the truth (1>0, 2>0, 2>1)
    false_negative = (
        matrix.at["1", "0"] if ("1" in matrix.index and "0" in matrix.columns) else 0
    ) + (
        matrix.at["2", "0"] if ("2" in matrix.index and "0" in matrix.columns) else 0
    ) + (
        matrix.at["2", "1"] if ("2" in matrix.index and "1" in matrix.columns) else 0
    )

    false_positive_rate = false_positive / total if total > 0 else None
    false_negative_rate = false_negative / total if total > 0 else None
    
    # --- VMISS statistics ---
    vmiss_mean = vmiss_df["VMISS_FREQ"].mean() if vmiss_df is not None else None
    vmiss_sd = vmiss_df["VMISS_FREQ"].std() if vmiss_df is not None else None
    vmiss_median = vmiss_df["VMISS_FREQ"].median() if vmiss_df is not None else None
    
    # --- SMISS statistics ---
    smiss_mean = smiss_df["SMISS_FREQ"].mean() if smiss_df is not None else None
    smiss_sd = smiss_df["SMISS_FREQ"].std() if smiss_df is not None else None
    smiss_median = smiss_df["SMISS_FREQ"].median() if smiss_df is not None else None

    # --- One-row summary ---
    summary_df = pd.DataFrame([{
        'DP': global_dp,
        'GQ': global_gq,
        'LAF': global_laf,
        'HAF': global_haf,
        'PLATFORM': platform,
        'TOTAL_GENOTYPE': int(total),
        'TOTAL_GENOTYPE(WITH_EMPTY)': int(total_with_empty),
        'TOTAL_GENOTYPE(WITH_ALL_NA)': int(matrix.values.sum()),
        'GENOTYPE_CONCORDANCE': concordance,
        'GENOTYPE_CONCORDANCE(WITH_EMPTY)': concordance_with_empty,
        'GENOTYPE_MISS_RATE': miss_rate,
        'FALSE_POSITIVE_RATE': false_positive_rate,
        'FALSE_NEGATIVE_RATE': false_negative_rate,
        'HET>HOMVAR_COUNT': int(het_to_homvar),
        'HOMREF>HOMVAR_COUNT': int(homref_to_homvar),
        'HOMREF>HET_COUNT': int(homref_to_het),
        'HET>HOMREF_COUNT': int(het_to_homref),
        'HOMVAR>HOMREF_COUNT': int(homvar_to_homref),
        'HOMVAR>HET_COUNT': int(homvar_to_het),
        'HOMREF>HOMREF_COUNT': int(homref_to_homref),
        'HET>HET_COUNT': int(het_to_het),
        'HOMVAR>HOMVAR_COUNT': int(homvar_to_homvar),
        'HET>HOMVAR_RATE': het_to_homvar_rate,
        'HOMREF>HOMVAR_RATE': homref_to_homvar_rate,
        'HOMREF>HET_RATE': homref_to_het_rate,
        'HET>HOMREF_RATE': het_to_homref_rate,
        'HOMVAR>HOMREF_RATE': homvar_to_homref_rate,
        'HOMVAR>HET_RATE': homvar_to_het_rate,
        'HOMREF>HOMREF_RATE': homref_to_homref_rate,
        'HET>HET_RATE': het_to_het_rate,
        'HOMVAR>HOMVAR_RATE': homvar_to_homvar_rate,
        'VMISS_FREQ_MEAN': vmiss_mean,
        'VMISS_FREQ_SD': vmiss_sd,
        'VMISS_FREQ_MEDIAN': vmiss_median,
        'SMISS_FREQ_MEAN': smiss_mean,
        'SMISS_FREQ_SD': smiss_sd,
        'SMISS_FREQ_MEDIAN': smiss_median
    }])

    return summary_df

# --- Per sample group ---
# For each sample group (ALL, 15X, 30X): take its subset, compute the summary and store it in the output
# dictionary.
platforms = ['ALL', '15X', '30X']
results_dict = {}

for platform in platforms:
    vmiss_sub = vmiss_df[vmiss_df["PLATFORM"] == platform]
    smiss_sub = smiss_df[smiss_df["PLATFORM"] == platform]
    confusion_sub = confusion_df[confusion_df["PLATFORM"] == platform]
    matrix = confusion_to_matrix(confusion_sub)
    summary = summarize_genotype_concordance(matrix, vmiss_sub, smiss_sub, global_dp, global_gq, global_laf, global_haf, platform)

    key = frozenset([f'DP{global_dp}', f'GQ{global_gq}', f'LAF{global_laf}', f'HAF{global_haf}', platform])
    results_dict[key] = (vmiss_sub, smiss_sub, confusion_sub, summary)

# --- Write ---
# Save the summary as a .summary.pkl file named after the thresholds.
summary_path = f'{chr}._DP{global_dp}_GQ{global_gq}_LAF{global_laf}_HAF{global_haf}_.summary.pkl'
with open(summary_path, 'wb') as f:
    pickle.dump(results_dict, f)
