/*
 * ParetoGT: Pareto-optimal genotype QC thresholds for WGS, calibrated against array genotypes.
 *
 * For every combination of a DP threshold, a GQ threshold and a heterozygote allele-fraction window
 * [LAF, HAF], WGS genotypes failing the thresholds are set to missing and compared with the array genotypes
 * of the same samples at the shared variants. The settings on the Pareto front of concordance against miss
 * rate are then ranked by a weighted Kneedle.
 *
 *   PREPARE_MATRICES     per chromosome: shared samples and variants; GT/DP/GQ/AF matrices      -> 01
 *                        (skipped with --matrixDir: the matrices of an earlier run)
 *   EVALUATE             per chromosome x batch of combinations: confusion matrix, VMISS, SMISS -> 02
 *                        (and the per-chromosome metrics, --chromSummaries)                     -> 03
 *   MERGE_AND_SUMMARIZE  per batch: genome-wide merge and metrics of each combination           -> 04, 05
 *   EXPORT_TABLES        the metrics of all combinations, one table per sample group           -> summary
 *   SELECT_SETTING       Pareto front and the knee point of each weight                        -> selection
 *   FIGURES              threshold scan, Pareto front, allele-fraction window, depth groups    -> figures
 *
 * Parameters and profiles: nextflow.config. Methods: docs/METHODS.md.
 */

// "1-3,5" -> ['1', '2', '3', '5']; values stay strings
def expandList(spec) {
    def out = []
    spec.toString().tokenize(',').each { token ->
        def t = token.trim()
        if (t ==~ /\d+-\d+/) {
            def ends = t.tokenize('-')
            (ends[0].toInteger()..ends[1].toInteger()).each { v -> out << v.toString() }
        } else if (t) {
            out << t
        }
    }
    return out
}

// Published per chromosome, as <chr>/<file>: the file names start with "<chr>."
def byChromosome(fn) {
    return "${fn.tokenize('.')[0]}/${fn}"
}

// Thresholds written as the Python scripts write them in file names: 8, 20, 0.2, 1.0
def asInt(v) {
    return v.toString().toInteger().toString()
}

def asFloat(v) {
    return v.toString().toDouble().toString()
}

process PREPARE_MATRICES {
    tag "${chr}"
    label 'prepare'
    publishDir "${params.outDir}/01.prepare_format_matrix", mode: params.publishMode, saveAs: { fn -> byChromosome(fn) }

    input:
    val chr
    path extract_py

    output:
    tuple val(chr), path("${chr}.shared_samples.txt"), path("${chr}.shared_variant_ids.txt"),
          path("${chr}.wgs.GT.tsv.gz"), path("${chr}.wgs.DP.tsv.gz"), path("${chr}.wgs.GQ.tsv.gz"),
          path("${chr}.wgs.AF.tsv.gz"), path("${chr}.array.GT.tsv.gz")

    script:
    def callVcf = "${params.wgsDir}/${chr}${params.callVcfSuffix}"
    """
    export PATH=\$(dirname ${params.bcftools}):\$PATH
    ${params.python} ${extract_py} \\
        --chr ${chr} \\
        --call_prefix wgs \\
        --true_prefix array \\
        --call_vcf ${callVcf} \\
        --true_vcf ${params.trueVcf} \\
        --num_threads ${params.prepareThreads} \\
        --call_format_fields GT DP GQ AF
    """
}

process EVALUATE {
    tag "${chr} batch ${batch}"
    label 'evaluate'
    publishDir "${params.outDir}/02.concordance_vmiss_calculator", mode: params.publishMode, pattern: 'chr*._*_.pkl', saveAs: { fn -> byChromosome(fn) }
    publishDir "${params.outDir}/03.concordance_vmiss_summary", mode: params.publishMode, pattern: 'chr*._*_.summary.pkl', saveAs: { fn -> byChromosome(fn) }

    input:
    tuple val(chr), path(samples), path(variants), path(gt), path(dp), path(gq), path(af), path(truth), val(batch), val(combos)
    path sample_info
    path evaluate_py
    path summary_py

    output:
    tuple val(batch), path("${chr}._*_.pkl"), emit: pkl
    path "${chr}._*_.summary.pkl", optional: true, emit: summary

    script:
    def comboArgs = combos.collect { c -> c.join(',') }.join(' ')
    def summaries = params.chromSummaries ? combos.collect { c ->
        "${params.python} ${summary_py} --chr ${chr} " +
        "--global_dp ${c[0]} --global_gq ${c[1]} --global_laf ${c[2]} --global_haf ${c[3]} " +
        "--raw_pkl ${chr}._DP${c[0]}_GQ${c[1]}_LAF${c[2]}_HAF${c[3]}_.pkl"
    }.join('\n') : ''
    """
    ${params.python} ${evaluate_py} \\
        --chr ${chr} \\
        --combos ${comboArgs} \\
        --batch_size ${params.batchSize} \\
        --num_threads ${params.evaluateThreads} \\
        --call_gt ${gt} \\
        --call_dp ${dp} \\
        --call_gq ${gq} \\
        --call_af ${af} \\
        --true_gt ${truth} \\
        --samples_txt ${samples} \\
        --variants_txt ${variants} \\
        --sample_info_csv ${sample_info}
    ${summaries}
    """
}

process MERGE_AND_SUMMARIZE {
    tag "batch ${batch}"
    label 'merge'
    publishDir "${params.outDir}/04.merge_concordance_vmiss", mode: params.publishMode, pattern: 'merged.*_.pkl'
    publishDir "${params.outDir}/05.merge_concordance_vmiss_summary", mode: params.publishMode, pattern: 'merged.*.summary.pkl'

    input:
    tuple val(batch), val(combos), path(chrom_pkls)
    path merge_py
    path summary_py

    output:
    path 'merged.*_.pkl', emit: merged
    path 'merged.*.summary.pkl', emit: summary

    script:
    def steps = combos.collect { c ->
        def name = "_DP${c[0]}_GQ${c[1]}_LAF${c[2]}_HAF${c[3]}_"
        "${params.python} ${merge_py} " +
        "--dp ${c[0]} --gq ${c[1]} --laf ${c[2]} --haf ${c[3]} " +
        "--input_chr_pkl \$(ls chr*.${name}.pkl | sort -V)\n" +
        "${params.python} ${summary_py} --chr merged " +
        "--global_dp ${c[0]} --global_gq ${c[1]} --global_laf ${c[2]} --global_haf ${c[3]} " +
        "--raw_pkl merged.${name}.pkl"
    }.join('\n')
    """
    ${steps}
    """
}

process EXPORT_TABLES {
    label 'export'
    publishDir "${params.outDir}/summary", mode: 'copy'

    input:
    path summaries
    val n_expected
    path export_py

    output:
    path 'merged_summary_*.csv'

    script:
    """
    ${params.python} ${export_py} \\
        --summary_dir . \\
        --out_dir . \\
        --expect ${n_expected} \\
        --workers ${params.exportThreads}
    """
}

process SELECT_SETTING {
    label 'small'
    publishDir "${params.outDir}/selection", mode: 'copy'

    input:
    path tables
    path select_py

    output:
    path '*.tsv'

    script:
    """
    ${params.python} ${select_py} --tables_dir . --out_dir . --weights '${params.weights}'
    """
}

process FIGURES {
    label 'small'
    publishDir "${params.outDir}/figures", mode: 'copy'

    input:
    path tables
    path selection
    path figures_py

    output:
    path '*.{png,pdf}'

    script:
    """
    export MPLCONFIGDIR=\$PWD/.mpl
    ${params.python} ${figures_py} \\
        --tables_dir . \\
        --knee_table knee_by_weight.tsv \\
        --out_dir . \\
        --adopt_weight '${params.adoptWeight}'
    """
}

workflow {
    def chroms = expandList(params.chromosomes).collect { c -> c.startsWith('chr') ? c : "chr${c}".toString() }
    def grid = []
    expandList(params.dp).each { dp ->
        expandList(params.gq).each { gq ->
            expandList(params.laf).each { laf ->
                expandList(params.haf).each { haf ->
                    grid << [asInt(dp), asInt(gq), asFloat(laf), asFloat(haf)]
                }
            }
        }
    }
    def batches = grid.collate(params.combosPerTask as int)
    def indexed = (0..<batches.size()).collect { i -> tuple(i, batches[i]) }
    log.info "ParetoGT: ${chroms.size()} chromosomes x ${grid.size()} threshold combinations, in ${batches.size()} batches"

    sample_info = file(params.sampleInfo, checkIfExists: true)
    summary_py = file("${projectDir}/scripts/genotype_concordance_vmiss_summary.py", checkIfExists: true)

    def prepared
    if (params.matrixDir) {
        // matrices of an earlier run (its 01.prepare_format_matrix): the VCFs are not read again
        prepared = channel.fromList(chroms).map { c ->
            def m = "${params.matrixDir}/${c}/${c}"
            tuple(c, file("${m}.shared_samples.txt", checkIfExists: true),
                  file("${m}.shared_variant_ids.txt", checkIfExists: true),
                  file("${m}.wgs.GT.tsv.gz", checkIfExists: true), file("${m}.wgs.DP.tsv.gz", checkIfExists: true),
                  file("${m}.wgs.GQ.tsv.gz", checkIfExists: true), file("${m}.wgs.AF.tsv.gz", checkIfExists: true),
                  file("${m}.array.GT.tsv.gz", checkIfExists: true))
        }
    } else {
        if (!params.wgsDir || !params.trueVcf) {
            error "Give the inputs: --wgsDir and --trueVcf (or --matrixDir, the matrices of an earlier run)"
        }
        PREPARE_MATRICES(channel.fromList(chroms), file("${projectDir}/scripts/extract_vcf_format.py", checkIfExists: true))
        prepared = PREPARE_MATRICES.out
    }
    EVALUATE(prepared.combine(channel.fromList(indexed)), sample_info,
             file("${projectDir}/scripts/evaluate_genotype_concordance_and_vmiss.py", checkIfExists: true), summary_py)

    per_batch = EVALUATE.out.pkl
        .groupTuple(size: chroms.size())
        .map { batch, pkls -> tuple(batch, batches[batch], pkls.flatten()) }
    MERGE_AND_SUMMARIZE(per_batch, file("${projectDir}/scripts/merge_concordance_results.py", checkIfExists: true), summary_py)

    EXPORT_TABLES(MERGE_AND_SUMMARIZE.out.summary.flatten().collect(), grid.size(),
                  file("${projectDir}/scripts/export_summary_tables.py", checkIfExists: true))
    SELECT_SETTING(EXPORT_TABLES.out, file("${projectDir}/scripts/select_setting.py", checkIfExists: true))
    FIGURES(EXPORT_TABLES.out, SELECT_SETTING.out, file("${projectDir}/scripts/plot_figures.py", checkIfExists: true))
}
