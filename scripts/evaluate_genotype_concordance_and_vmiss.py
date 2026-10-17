"""
Script: evaluate_genotype_concordance_and_vmiss.py
Purpose: evaluate the concordance of WGS genotypes with array genotypes, and compute the missing rates per
sample and per variant (SMISS and VMISS), under one or more combinations of genotype QC thresholds.

Description:
    Reads the WGS genotype (GT), read depth (DP), genotype quality (GQ) and allele fraction (AF) matrices and
    the array genotypes of one chromosome, once. For each threshold combination, sets a WGS genotype to
    missing when it fails the thresholds (DP, GQ; for a heterozygous call, AF outside [LAF, HAF]), then builds
    the WGS x array confusion matrix and counts the missing genotypes per variant (VMISS) and per sample
    (SMISS), in all samples and in the 15x and 30x samples.

    The variants are processed in blocks of --batch_size, in parallel, exactly as in the first version of
    this script (which re-read the matrices for every block); only the reading is shared now, so the output
    of a combination is identical.

Inputs:
    - {chr}.wgs.GT.tsv.gz: WGS genotypes
    - {chr}.wgs.DP.tsv.gz: WGS read depth
    - {chr}.wgs.GQ.tsv.gz: WGS genotype quality
    - {chr}.wgs.AF.tsv.gz: WGS allele fraction
    - {chr}.array.GT.tsv.gz: array genotypes (truth)
    - {chr}.shared_samples.txt: shared samples
    - {chr}.shared_variant_ids.txt: shared variants
    - sample information CSV (ID, Target_DP): sequencing depth of each sample (15x or 30x)

Thresholds:
    --combos DP,GQ,LAF,HAF [DP,GQ,LAF,HAF ...]   one or more combinations; or, for one combination,
    --global_dp, --global_gq, --global_low_af, --global_high_af

Output, one file per combination:
    {chr}._DP{dp}_GQ{gq}_LAF{laf}_HAF{haf}_.pkl: a dictionary keyed by the combination and sample group
    (ALL / 15X / 30X), holding (vmiss_df, smiss_df, confusion_df)
        * vmiss_df: missing count and missing rate of each variant
        * smiss_df: missing count and missing rate of each sample
        * confusion_df: WGS x array confusion matrix
"""
import argparse
import multiprocessing as mp
import pickle
from functools import partial

import pandas as pd

# === Arguments ===
parser = argparse.ArgumentParser(description="Evaluate genotype concordance and VMISS for WGS vs array data.")
parser.add_argument("--chr", type=str, required=True, help="Chromosome to process (e.g., chr22)")
parser.add_argument("--combos", type=str, nargs="+", default=None,
                    help="Threshold combinations DP,GQ,LAF,HAF (overrides the --global_* thresholds)")
parser.add_argument("--global_dp", type=int, default=8, help="Global DP threshold for filtering")
parser.add_argument("--global_gq", type=int, default=20, help="Global GQ threshold for filtering")
parser.add_argument("--global_low_af", type=float, default=0.25, help="Low allele fraction threshold for filtering")
parser.add_argument("--global_high_af", type=float, default=0.75, help="High allele fraction threshold for filtering")
parser.add_argument("--batch_size", type=int, default=1000, help="Batch size for processing variants")
parser.add_argument("--num_threads", type=int, default=4, help="Number of threads for multiprocessing")
parser.add_argument("--call_gt", type=str, required=True, help="Path to WGS genotype file (GT)")
parser.add_argument("--call_dp", type=str, required=True, help="Path to WGS depth file (DP)")
parser.add_argument("--call_gq", type=str, required=True, help="Path to WGS genotype quality file (GQ)")
parser.add_argument("--call_af", type=str, required=True, help="Path to WGS allele fraction file (AF)")
parser.add_argument("--true_gt", type=str, required=True, help="Path to true genotype file (GT) for array data")
parser.add_argument("--samples_txt", type=str, required=True, help="Path to shared samples text file")
parser.add_argument("--variants_txt", type=str, required=True, help="Path to shared variant IDs text file")
parser.add_argument("--sample_info_csv", type=str, required=True, help="Path to sample info CSV file with DP information")
args = parser.parse_args()

chr = args.chr
batch_size = args.batch_size  # Batch size for processing variants
num_threads = args.num_threads  # Number of threads for multiprocessing

if args.combos:
    combos = []
    for c in args.combos:
        dp, gq, laf, haf = c.split(",")
        combos.append((int(dp), int(gq), float(laf), float(haf)))
else:
    combos = [(args.global_dp, args.global_gq, args.global_low_af, args.global_high_af)]

# === Read the sample and variant IDs ===
with open(args.samples_txt) as f:
    sample_ids = [line.strip() for line in f if line.strip()]
with open(args.variants_txt) as f:
    variant_ids = [line.strip() for line in f if line.strip()]

# === Read the sequencing depth of each sample ===
sample_info = pd.read_csv(args.sample_info_csv, usecols=["ID", "Target_DP"])

# === Sample IDs at 15x and 30x, restricted to sample_ids ===
sample_ids_set = set(sample_ids)

sample_ids_15x = sample_info.loc[sample_info["Target_DP"] == "15x", "ID"]
sample_ids_30x = sample_info.loc[sample_info["Target_DP"] == "30x", "ID"]

sample_ids_15x = list(sample_ids_set & set(sample_ids_15x))
sample_ids_30x = list(sample_ids_set & set(sample_ids_30x))

# === The matrices, read once (samples x variants, in file order) ===
CALL_GT = pd.read_csv(args.call_gt, sep="\t", compression="gzip", index_col=0)
CALL_DP = pd.read_csv(args.call_dp, sep="\t", compression="gzip", index_col=0)
CALL_GQ = pd.read_csv(args.call_gq, sep="\t", compression="gzip", index_col=0)
CALL_AF = pd.read_csv(args.call_af, sep="\t", compression="gzip", index_col=0)
TRUE_GT = pd.read_csv(args.true_gt, sep="\t", compression="gzip", index_col=0)

# Thresholds of the combination being evaluated; set before each pool is forked
global_dp = global_gq = global_low_af = global_high_af = None


# === Top-level function, run in worker processes ===
# - takes the GT/DP/GQ/AF matrices of a sample subset and a variant block
# - sets genotypes failing the DP/GQ/AF thresholds to missing
# - from the filtered genotypes:
#     * builds the confusion matrix (confusion_df)
#     * computes the per-variant missing counts (vmiss_df)
#     * computes the per-sample missing counts (smiss_df)
def process_variant_block(cols, sample_ids_subset):
    gt = CALL_GT.loc[sample_ids_subset, cols].astype("Int64")
    dp = CALL_DP.loc[sample_ids_subset, cols]
    gq = CALL_GQ.loc[sample_ids_subset, cols]
    af = CALL_AF.loc[sample_ids_subset, cols]
    tru = TRUE_GT.loc[sample_ids_subset, cols].astype("Int64")

    assert list(gt.columns) == list(tru.columns), "Variant ID columns not aligned"
    assert list(gt.index) == list(tru.index), "Sample ID rows not aligned"

    qc_gt = gt.copy()
    mask = (dp < global_dp) | (gq < global_gq)
    mask |= ((gt == 1) & (af < global_low_af))
    mask |= ((gt == 1) & (af > global_high_af))
    qc_gt[mask] = pd.NA

    merged = pd.DataFrame({
        "CALL_GENOTYPE": qc_gt.values.ravel(),
        "TRUE_GENOTYPE": tru.values.ravel()
    })
    counts = merged.groupby(["CALL_GENOTYPE", "TRUE_GENOTYPE"], dropna=False).size().reset_index(name="COUNT")

    vmiss_count = qc_gt.isna().sum()
    total_sample = qc_gt.shape[0]
    vmiss = pd.DataFrame({
        "VARIANT_ID": vmiss_count.index,
        "VMISS_COUNT": vmiss_count.values,
        "TOTAL_SAMPLE": total_sample
    })

    smiss_count = qc_gt.isna().sum(axis=1)
    total_variant = qc_gt.shape[1]
    smiss = pd.DataFrame({
        "SAMPLE_ID": smiss_count.index,
        "SMISS_COUNT": smiss_count.values,
        "TOTAL_VARIANT": total_variant
    })

    return counts, vmiss, smiss


# === Driver ===
# - runs process_variant_block for each sample group (ALL / 15X / 30X) and aggregates the blocks
# - sorts the variants by position and aggregates the missing counts per sample
# - stores each group in results_dict, keyed by the threshold combination and the group name
def run_and_save(results_dict, name, sample_ids):
    block_func = partial(process_variant_block, sample_ids_subset=sample_ids)
    with mp.get_context("fork").Pool(processes=num_threads) as pool:
        results = pool.map(block_func, blocks)

    # === Combine the blocks ===
    confusion_tables, vmiss_tables, smiss_tables = zip(*results)
    confusion_df = pd.concat(confusion_tables, ignore_index=True)
    confusion_df = confusion_df.groupby(["CALL_GENOTYPE", "TRUE_GENOTYPE"], dropna=False).sum().reset_index()

    def format_value(x):
        return "NA" if pd.isna(x) else str(int(x))
    confusion_df["CALL_GENOTYPE"] = confusion_df["CALL_GENOTYPE"].map(format_value)
    confusion_df["TRUE_GENOTYPE"] = confusion_df["TRUE_GENOTYPE"].map(format_value)

    vmiss_df = pd.concat(vmiss_tables, ignore_index=True)
    vmiss_df = vmiss_df.groupby("VARIANT_ID", as_index=False).agg({
        "VMISS_COUNT": "sum",
        "TOTAL_SAMPLE": "first"
    })
    vmiss_df["VMISS_FREQ"] = (vmiss_df["VMISS_COUNT"] / vmiss_df["TOTAL_SAMPLE"]).round(4)
    # === Sort the variants by position ===
    vmiss_df["__POS_INT"] = vmiss_df["VARIANT_ID"].str.split(":").str[1].astype(int)
    vmiss_df.sort_values("__POS_INT", inplace=True)
    vmiss_df.drop(columns="__POS_INT", inplace=True)
    vmiss_df.reset_index(drop=True, inplace=True)

    smiss_df = pd.concat(smiss_tables, ignore_index=True)
    smiss_df = smiss_df.groupby("SAMPLE_ID", as_index=False).agg({
        "SMISS_COUNT": "sum",
        "TOTAL_VARIANT": "sum"
    })
    smiss_df["SMISS_FREQ"] = (smiss_df["SMISS_COUNT"] / smiss_df["TOTAL_VARIANT"]).round(4)

    # === Key ===
    param_key = frozenset([
        f'DP{global_dp}',
        f'GQ{global_gq}',
        f'LAF{global_low_af}',
        f'HAF{global_high_af}',
        name.upper()  # ALL / 15X / 30X
    ])

    # === Store in the result dictionary ===
    results_dict[param_key] = (vmiss_df, smiss_df, confusion_df)
    print(f"[✓] Results stored in results_dict under key → {param_key}")


# === Variant blocks ===
blocks = [variant_ids[i:i + batch_size] for i in range(0, len(variant_ids), batch_size)]

# === Each combination: all, 15x and 30x, then one .pkl ===
for global_dp, global_gq, global_low_af, global_high_af in combos:
    results_dict = {}
    run_and_save(results_dict, "all", sample_ids)
    run_and_save(results_dict, "15x", sample_ids_15x)
    run_and_save(results_dict, "30x", sample_ids_30x)
    with open(f'{chr}._DP{global_dp}_GQ{global_gq}_LAF{global_low_af}_HAF{global_high_af}_.pkl', 'wb') as f:
        pickle.dump(results_dict, f)
