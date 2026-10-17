"""
extract_vcf_format.py
Extracts FORMAT-field matrices for the variants and samples shared by two VCF files (for example, WGS and
array) on one chromosome, as input to a genotype-concordance evaluation.
What it does:
1. Reads, in parallel, the variant IDs of both VCF files on the chromosome and keeps the shared ones.
2. Reads the sample IDs of both VCF files and keeps the shared ones.
3. Writes the shared variant IDs and sample IDs to text files.
4. Extracts, chunk by chunk with bcftools, the requested FORMAT fields (e.g. GT, DP, GQ, AF) of the shared
   variants and samples, and writes each as a compressed sample x variant matrix (tsv.gz).
Arguments:
- --call_prefix: prefix of the output files of the call VCF (e.g. wgs).
- --true_prefix: prefix of the output files of the truth VCF (e.g. array).
- --call_vcf: path of the call VCF (required).
- --true_vcf: path of the truth VCF (required).
- --chr: chromosome to process (required).
- --max_variants: maximum number of variants to process (for debugging; default: all).
- --num_threads: bcftools threads (default 4).
- --call_format_fields: FORMAT fields to extract from the call VCF (default ["GT", "DP", "GQ", "AF"]).
Outputs:
- {chr}.shared_variant_ids.txt: the shared variant IDs.
- {chr}.shared_samples.txt: the shared sample IDs.
- {chr}.{call_prefix}.{field}.tsv.gz: sample x variant matrix of each FORMAT field of the call VCF.
- {chr}.{true_prefix}.GT.tsv.gz: sample x variant matrix of GT in the truth VCF.
Requirements:
- bcftools
- numpy
- Python standard library (concurrent.futures, subprocess, argparse, tempfile, os, gc, shutil)
Example:
python extract_vcf_format.py --call_vcf call.vcf.gz --true_vcf truth.vcf.gz --chr chr1 --call_prefix wgs --true_prefix array --num_threads 8 --call_format_fields GT DP GQ AF
"""
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import subprocess
from typing import List
import subprocess
import argparse
import tempfile
import os
import gc
import shutil
import gzip

# %%
# Arguments
parser = argparse.ArgumentParser(description="Extract VCF format fields for shared variants and samples.")
parser.add_argument("--call_prefix", type=str, default="wgs", help="Prefix for call VCF files (default: 'wgs').")
parser.add_argument("--true_prefix", type=str, default="array", help="Prefix for true VCF files (default: 'array').")
parser.add_argument("--call_vcf", type=str, required=True, help="Path to the call VCF file.")
parser.add_argument("--true_vcf", type=str, required=True, help="Path to the true VCF file.")
parser.add_argument("--chr", type=str, required=True, help="Chromosome to process.")
parser.add_argument("--max_variants", type=int, default=None, help="Maximum number of variants to process (default: None, process all).")
parser.add_argument(
    "--num_threads", 
    type=int, 
    default=4, 
    help="Number of threads to use for bcftools (default: 4)."
)
parser.add_argument(
    "--call_format_fields", 
    type=str, 
    nargs='+', 
    default=["GT", "DP", "GQ", "AF"],
    help="List of FORMAT fields to extract (e.g. --call_format_fields GT DP GQ AF)"
)
args = parser.parse_args()

# %%
call_prefix = args.call_prefix
true_prefix = args.true_prefix

call_vcf_path = args.call_vcf
true_vcf_path = args.true_vcf

chr = args.chr

# %%

from typing import Optional

def extract_variant_ids(vcf_path: str, chr: str, max_variants: Optional[int] = None) -> set:
    """Return the set of ID fields of the VCF on chr."""
    query_cmd = [
        "bcftools", "query",
        "-r", chr,
        "-f", "%ID\n",
        vcf_path
    ]
    proc = subprocess.Popen(query_cmd, stdout=subprocess.PIPE, text=True, bufsize=8192)

    variant_ids = set()
    if proc.stdout is not None:
        for i, line in enumerate(proc.stdout):
            if max_variants is not None and i >= max_variants:
                break
            line = line.rstrip('\n')
            if line != ".":
                variant_ids.add(line)
        proc.stdout.close()
    proc.wait()
    return variant_ids


def extract_sample_ids(vcf_path: str) -> list:
    """Return the list of sample IDs of the VCF."""
    header = subprocess.check_output(["bcftools", "view", "-h", vcf_path], text=True)
    for line in header.strip().splitlines():
        if line.startswith("#CHROM"):
            return line.strip().split("\t")[9:]
    return []


# === Settings ===

max_variants = args.max_variants

# === Read the IDs in parallel ===
from concurrent.futures import wait

with ThreadPoolExecutor(max_workers=2) as executor:
    futures = {
        "true": executor.submit(extract_variant_ids, true_vcf_path, chr, max_variants),
        "call": executor.submit(extract_variant_ids, call_vcf_path, chr, max_variants)
    }
    done, not_done = wait(futures.values())
    for fut in done:
        exc = fut.exception()
        if exc is not None:
            raise exc
    true_ids = futures["true"].result()
    call_ids = futures["call"].result()

# === Shared variant IDs and shared samples ===
shared_variant_ids = true_ids & call_ids
shared_variant_ids = sorted(shared_variant_ids, key=lambda x: int(x.split(":")[1]))
true_samples = extract_sample_ids(true_vcf_path)
call_samples = extract_sample_ids(call_vcf_path)
shared_samples = sorted(set(true_samples) & set(call_samples))  # deterministic order


# === Report the counts ===
print(f"Shared samples: {len(shared_samples)}")
print(f"Shared variant IDs on {chr}: {len(shared_variant_ids)}")

# === Free memory explicitly ===
del true_ids, call_ids, true_samples, call_samples
gc.collect()

# %%
# === Write the shared variant IDs and sample IDs ===

with open(f"{chr}.shared_variant_ids.txt", "w") as f:
    for vid in shared_variant_ids:
        f.write(vid + "\n")

with open(f"{chr}.shared_samples.txt", "w") as f:
    for sid in shared_samples:
        f.write(sid + "\n")

# %%
def extract_format_matrices_to_disk(
    vcf_path: str,
    sample_ids: List[str],
    focus_loci: List[str],
    format_fields_to_extract: List[str] = ["GT"],
    use_gt_map: bool = True,
    num_threads: int = 4
) -> dict:
    """
    Extract the FORMAT-field matrices of the given samples and variants from a VCF file.

    Args:
        vcf_path (str): path of the input VCF.
        sample_ids (List[str]): sample IDs to extract.
        focus_loci (List[str]): variant IDs to extract, as 'chr:pos:ref:alt'.
        format_fields_to_extract (List[str], optional): FORMAT fields to extract; default ["GT"].
        use_gt_map (bool, optional): encode GT as numbers; default True.
        num_threads (int, optional): bcftools threads; default 4.

    Returns:
        dict: {field_name: {sample_id: [values in locus order]}}
    """
    tmp_dir = tempfile.mkdtemp(prefix="matrix_tmp_")

    # VCF column indices
    VCF_COLUMNS = {
        "CHROM": 0,
        "POS": 1,
        "ID": 2,
        "REF": 3,
        "ALT": 4,
        "QUAL": 5,
        "FILTER": 6,
        "INFO": 7,
        "FORMAT": 8,
        "SAMPLES_START": 9
    }

    # GT encoded as numbers
    gt_map = {
        "0/0": 0, "0|0": 0,
        "0/1": 1, "1/0": 1, "0|1": 1, "1|0": 1,
        "1/1": 2, "1|1": 2,
        "./.": np.nan, ".|.": np.nan
    }

    # Region list for the query and a lookup of the loci
    loci_dict = {}
    for locus in focus_loci:
        chrom, pos, ref, alt = locus.split(":")
        key = (chrom, int(pos), ref, alt)
        loci_dict[key] = locus
    
    region_file_path = os.path.join(tmp_dir, "regions.txt")
    with open(region_file_path, "w") as f:
        for locus in focus_loci:
            chrom, pos, ref, alt = locus.split(":")
            f.write(f"{chrom}\t{pos}\t{pos}\n")
    
    # Sample names of the VCF
    header = subprocess.check_output(["bcftools", "view", "-h", vcf_path], text=True)
    vcf_samples = None
    for line in header.strip().splitlines():
        if line.startswith("#CHROM"):
            vcf_samples = line.strip().split('\t')[VCF_COLUMNS["SAMPLES_START"]:]
            break
    if vcf_samples is None:
        raise RuntimeError("Could not find #CHROM header line in VCF to extract sample IDs.")

    # Result dictionary
    result = {field: {s: [] for s in sample_ids} for field in format_fields_to_extract}

    # bcftools command
    cmd = [
        "bcftools", "view",
        "--threads", str(num_threads),
        "-R", region_file_path,
        vcf_path
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)

    if proc.stdout is not None:
        for line in proc.stdout:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[VCF_COLUMNS["CHROM"]]
            pos = int(fields[VCF_COLUMNS["POS"]])
            ref = fields[VCF_COLUMNS["REF"]]
            alt = fields[VCF_COLUMNS["ALT"]].split(',')[0]
            key = (chrom, pos, ref, alt)
            if key not in loci_dict:
                continue
            locus_id = loci_dict[key]

            format_keys = fields[VCF_COLUMNS["FORMAT"]].split(":")
            format_index = {k: i for i, k in enumerate(format_keys)}
            sample_fields = fields[VCF_COLUMNS["SAMPLES_START"]:]

            for field in format_fields_to_extract:
                values = []
                for sample_str in sample_fields:
                    parts = sample_str.split(":")
                    raw_value = parts[format_index[field]] if field in format_index and len(parts) > format_index[field] else None
                    if field == "GT":
                        key = raw_value if raw_value is not None else "."
                        val = gt_map.get(key, np.nan) if use_gt_map else raw_value
                    elif field == "AD":
                        try:
                            val = tuple(map(int, raw_value.split(","))) if raw_value and "," in raw_value else np.nan
                        except:
                            val = np.nan
                    else:
                        try:
                            val = float(raw_value) if (raw_value is not None and raw_value != ".") else np.nan
                        except:
                            val = np.nan
                    values.append(val)
                per_sample = dict(zip(vcf_samples, values))
                for s in sample_ids:
                    result[field][s].append(per_sample.get(s, np.nan))

        proc.stdout.close()
    proc.wait()

    shutil.rmtree(tmp_dir)

    return result

def write_matrix_chunk_v2(output_path, sample_ids, loci_chunk, data_dict, append=False):
    mode = "at" if append else "wt"
    with gzip.open(output_path, mode) as out:
        if not append:
            out.write("sample_id\t" + "\t".join(loci_chunk) + "\n")
        for sid in sample_ids:
            row = [sid]
            values = data_dict.get(sid, [])
            for v in values:
                if isinstance(v, tuple):
                    val_str = ",".join(map(str, v))
                elif v is np.nan or (isinstance(v, float) and np.isnan(v)):
                    val_str = ""
                else:
                    val_str = str(v)
                row.append(val_str)
            out.write("\t".join(row) + "\n")

def merge_matrix_chunks(output_path, sample_ids, all_loci, chunk_paths):
    with gzip.open(output_path, "wt") as out:
        out.write("sample_id\t" + "\t".join(all_loci) + "\n")  # header: all loci
        for sid in sample_ids:
            row = [sid]
            for tmp_path in chunk_paths:
                with gzip.open(tmp_path, "rt") as f:
                    header = f.readline()
                    found = False
                    for line in f:
                        if line.startswith(sid + "\t"):  # the row of this sample
                            row.extend(line.strip().split("\t")[1:])  # append its values
                            found = True
                            break
                    if not found:
                        row.extend([""] * (len(header.strip().split("\t")) - 1))
            out.write("\t".join(row) + "\n")  # one row per sample


# === Chunks processed in parallel, then merged ===
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_vcf_in_chunks_parallel(vcf_path, sample_ids, variant_ids, output_prefix, format_fields_to_extract, num_threads):
    chunk_size = 10000
    n = len(variant_ids)
    tmp_dir = tempfile.mkdtemp(prefix="matrix_chunk_")

    chunk_files_by_field = {field: {} for field in format_fields_to_extract}
    loci_chunks = {}

    def process_chunk(i, loci_chunk):
        data = extract_format_matrices_to_disk(
            vcf_path=vcf_path,
            sample_ids=sample_ids,
            focus_loci=loci_chunk,
            format_fields_to_extract=format_fields_to_extract,
            use_gt_map=True,
            num_threads=num_threads
        )
        paths = {}
        for field in format_fields_to_extract:
            tmp_path = os.path.join(tmp_dir, f"{field}.chunk{i}.tmp.gz")
            write_matrix_chunk_v2(tmp_path, sample_ids, loci_chunk, data[field], append=False)
            paths[field] = tmp_path
        return i, loci_chunk, paths

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for i in range(0, n, chunk_size):
            loci_chunk = variant_ids[i:i+chunk_size]
            loci_chunks[i] = loci_chunk
            futures.append(executor.submit(process_chunk, i, loci_chunk))
        for future in as_completed(futures):
            i, loci_chunk, paths = future.result()
            for field in format_fields_to_extract:
                chunk_files_by_field[field][i] = paths[field]

    def merge_chunks_for_field(field):
        output_path = f"{chr}.{output_prefix}.{field}.tsv.gz"
        sorted_indices = sorted(chunk_files_by_field[field].keys())
        chunk_paths = [chunk_files_by_field[field][i] for i in sorted_indices]
        for tmp_path in chunk_paths:
            if not os.path.exists(tmp_path):
                raise FileNotFoundError(f"Missing chunk file: {tmp_path}")
        all_loci = []
        for i in sorted_indices:
            all_loci.extend(loci_chunks[i])
        merge_matrix_chunks(output_path, sample_ids, all_loci, chunk_paths)

    merge_threads = min(num_threads, len(format_fields_to_extract))
    with ThreadPoolExecutor(max_workers=merge_threads) as merge_executor:
        merge_futures = [merge_executor.submit(merge_chunks_for_field, field) for field in format_fields_to_extract]
        for f in as_completed(merge_futures):
            f.result()

    shutil.rmtree(tmp_dir)

# Process call_vcf and true_vcf concurrently
with ThreadPoolExecutor(max_workers=2) as outer_executor:
    f_call = outer_executor.submit(
        process_vcf_in_chunks_parallel,
        call_vcf_path,
        shared_samples,
        shared_variant_ids,
        call_prefix,
        args.call_format_fields,
        args.num_threads
    )
    f_true = outer_executor.submit(
        process_vcf_in_chunks_parallel,
        true_vcf_path,
        shared_samples,
        shared_variant_ids,
        true_prefix,
        ["GT"],
        args.num_threads
    )
    f_call.result()
    f_true.result()