// =============================================================================
// wgs_gt_qc_calibration: calibration of WGS genotype QC thresholds against array genotypes
//
// For every combination of a DP threshold (1-30), a GQ threshold (10, 20, 30) and a heterozygote
// allele-fraction window [LAF, HAF] (LAF 0-0.3, HAF 0.7-1.0), WGS genotypes failing the thresholds are set
// to missing and compared with the array genotypes of the same samples at the shared variants:
//   01 FormatMatrixPrepare           shared samples and variants; GT/DP/GQ/AF matrices (per chromosome)
//   02 ConcordanceVmissCalculator    confusion matrix, VMISS and SMISS (per chromosome x combination)
//   03 ConcordanceVmissSummary       metrics per chromosome x combination
//   04 MergeConcordanceVmiss         genome-wide merge (per combination)
//   05 MergeConcordanceVmissSummary  genome-wide metrics (per combination): the input of the notebooks
// 22 chromosomes x 3,240 combinations. DSL1: run with Nextflow 22.10.x (conda env dsl1) through
// run_nextflow.sh. Selection of a setting (Pareto front, weighted Kneedle): notebooks/.
// =============================================================================

// ---- Inputs (site data; override on the command line, e.g. --wgsDir <dir>)
params.wgsDir        = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/wgs/results/03.af_gtnorm"    // WGS calls, one VCF per chromosome
params.callVcfSuffix = ".vqc.af.gtnorm.vcf.gz"                                                                        // <chr><suffix> in wgsDir
params.arrayDir      = "/LARGE0/gr10478/b37974/Pulmonary_Hypertension/cteph_agp3k.v5/array/results/14.final_vcf"  // array genotypes (truth)
params.trueVcf       = "${params.arrayDir}/cteph_agp3k.array.sqc.vqc.tommo.tommo.vcf.gz"
params.infoDir       = "${projectDir}/info"
params.sampleInfo    = "${params.infoDir}/jhrpv5.target_dp.csv"   // ID,Target_DP (15x / 30x) of every sample

// ---- This repository
params.scriptDir     = "${projectDir}/scripts"
params.outDir        = "${projectDir}/results"

Channel
    .from((1..22).collect { "chr${it}" })
    .set { chr_ch }

// chr_ch = Channel.from([1, 22])
//                 .map { it -> "chr${it}" }

process FormatMatrixPrepare {

    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${chr}"

    publishDir "${params.outDir}/01.prepare_format_matrix/${chr}", mode: 'symlink'

    input:
    val(chr) from chr_ch
    val(wgsDir) from params.wgsDir
    val(arrayDir) from params.arrayDir

    output:
    tuple val(chr),
          file("${chr}.shared_samples.txt"),
          file("${chr}.shared_variant_ids.txt"),
          file("${chr}.${call_prefix}.GT.tsv.gz"),
          file("${chr}.${call_prefix}.DP.tsv.gz"),
          file("${chr}.${call_prefix}.GQ.tsv.gz"),
          file("${chr}.${call_prefix}.AF.tsv.gz"),
          file("${chr}.${true_prefix}.GT.tsv.gz") \
          into full_format_ch

    script:
    call_prefix = "wgs"
    true_prefix = "array"
    call_vcf = "${wgsDir}/${chr}${params.callVcfSuffix}"
    call_vcf_tbi = "${call_vcf}.tbi"
    true_vcf = params.trueVcf
    true_vcf_tbi = "${true_vcf}.tbi"
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/extract_vcf_format.py \
        --chr ${chr} \
        --call_prefix ${call_prefix} \
        --true_prefix ${true_prefix} \
        --call_vcf ${call_vcf} \
        --true_vcf ${true_vcf} \
        --num_threads 8 \
        --call_format_fields GT DP GQ AF
    """
}

// full_format_ch.view()

dp_channel = Channel.from(1..30)
gq_channel = Channel.of(10, 20, 30)
laf_channel = Channel.of(0.0, 0.1, 0.15, 0.2, 0.25, 0.3)
haf_channel = Channel.of(1.0, 0.9, 0.85, 0.8, 0.75, 0.7)


dp_channel
    .combine(gq_channel)
    .combine(laf_channel)
    .combine(haf_channel)
    .set { gt_tuning_ch } 

full_format_ch
    .combine(gt_tuning_ch)
    .set { full_format_with_tuning_ch }


process ConcordanceVmissCalculator {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${chr}: DP${dp} GQ${gq} LAF${laf} HAF${haf}"

    publishDir "${params.outDir}/02.concordance_vmiss_calculator/${chr}", mode: 'symlink'

    input:
    tuple val(chr), 
          file(sample_ls), 
          file(variant_ls), 
          file(call_gt), 
          file(call_dp), 
          file(call_gq), 
          file(call_af), 
          file(true_gt), 
          val(dp),
          val(gq),
          val(laf),
          val(haf) \
          from full_format_with_tuning_ch
    val(infoDir) from params.infoDir

    output:
    tuple val(chr), 
          val(dp), 
          val(gq), 
          val(laf), 
          val(haf), 
          file("${chr}._DP${dp}_GQ${gq}_LAF${laf}_HAF${haf}_.pkl") \
          into concordance_vmiss_ch, concordance_vmiss_ch_2

    script:
    sample_info = params.sampleInfo
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/evaluate_genotype_concordance_and_vmiss.py \
        --chr ${chr} \
        --global_dp ${dp} \
        --global_gq ${gq} \
        --global_low_af ${laf} \
        --global_high_af ${haf} \
        --batch_size 1000 \
        --num_threads 4 \
        --call_gt ${call_gt} \
        --call_dp ${call_dp} \
        --call_gq ${call_gq} \
        --call_af ${call_af} \
        --true_gt ${true_gt} \
        --samples_txt ${sample_ls} \
        --variants_txt ${variant_ls} \
        --sample_info_csv ${sample_info} \
    """
}

process ConcordanceVmissSummary {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "${chr}"

    publishDir "${params.outDir}/03.concordance_vmiss_summary/${chr}", mode: 'symlink'

    input:
    tuple val(chr), 
          val(dp), 
          val(gq), 
          val(laf), 
          val(haf), 
          file(pkl_file) \
          from concordance_vmiss_ch

    output:
    tuple val(chr), 
          val(dp), 
          val(gq), 
          val(laf), 
          val(haf), 
          file("${chr}._DP${dp}_GQ${gq}_LAF${laf}_HAF${haf}_.summary.pkl") \
          into concordance_vmiss_summary_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/genotype_concordance_vmiss_summary.py \
        --chr ${chr} \
        --global_dp ${dp} \
        --global_gq ${gq} \
        --global_laf ${laf} \
        --global_haf ${haf} \
        --raw_pkl ${pkl_file} \
    """
}


// --- Group concordance_vmiss_ch_2 by dp, gq, laf, haf, and aggregate pkl files sorted by chr1~chr22 ---

concordance_vmiss_ch_2
    .map { chr, dp, gq, laf, haf, pkl_file ->
        tuple(dp, gq, laf, haf, chr, pkl_file)
    }
    .groupTuple(by: [0,1,2,3])
    .set { grouped_ch }

process MergeConcordanceVmiss {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "DP${dp} GQ${gq} LAF${laf} HAF${haf}"

    publishDir "${params.outDir}/04.merge_concordance_vmiss", mode: 'symlink'

    input:
    tuple val(dp), val(gq), val(laf), val(haf), file(chr), file(pkl_files) from grouped_ch

    output:
    tuple val('merged'), val(dp), val(gq), val(laf), val(haf), file("*.pkl") into merge_concordance_vmiss_ch

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/merge_concordance_results.py \
        --dp ${dp} \
        --gq ${gq} \
        --laf ${laf} \
        --haf ${haf} \
        --input_chr_pkl ${pkl_files.join(' ')}
    """
}

process MergeConcordanceVmissSummary {
    executor 'slurm'
    queue 'gr10478b'
    time '6d'

    tag "DP${dp} GQ${gq} LAF${laf} HAF${haf}"

    publishDir "${params.outDir}/05.merge_concordance_vmiss_summary", mode: 'symlink'

    input:
    tuple val(tag), val(dp), val(gq), val(laf), val(haf), file(pkl_file) from merge_concordance_vmiss_ch

    output:
    file("*.summary.pkl")

    script:
    """
    source activate cteph_geno_pro
    python ${params.scriptDir}/genotype_concordance_vmiss_summary.py \
        --chr ${tag} \
        --global_dp ${dp} \
        --global_gq ${gq} \
        --global_laf ${laf} \
        --global_haf ${haf} \
        --raw_pkl ${pkl_file} \
    """
}

