# Sample sheet

`params.sampleInfo` (default `info/jhrpv5.target_dp.csv`) assigns each sample to a sequencing-depth group. Comma-separated, with a header:

| Column | Content |
|---|---|
| `ID` | sample ID, as in the VCF files |
| `Target_DP` | `15x` or `30x` |

The metrics are computed in all shared samples (ALL) and in each group (15X, 30X).

The sheet holds sample IDs, so it is kept out of version control (`.gitignore`).
