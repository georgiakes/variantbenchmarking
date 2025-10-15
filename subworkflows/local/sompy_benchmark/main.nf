//
// SOMATIC: SUBWORKFLOW FOR SOMPY BENCHMARKING
//

include { HAPPY_SOMPY          } from '../../../modules/nf-core/happy/sompy'
include { SPLIT_SOMPY_FEATURES } from '../../../modules/local/custom/split_sompy_features'
include { FILTER_VARIANTS      } from '../../../modules/local/custom/filter_variants/'

workflow SOMPY_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index,  regionsbed, targets_bed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    falsepositive_bed  // reference channel [val(meta), bed]
    ambiguous_beds     // reference channel [val(meta), bed]

    main:

    versions            = Channel.empty()
    tagged_variants_csv = Channel.empty()

    // apply sompy for small somatic variant benchmarking
    HAPPY_SOMPY(
        input_ch.map{meta, test, test_index, truth, truth_index, regions, target -> [meta, test, truth, regions, target]},
        fasta,
        fai,
        falsepositive_bed,
        ambiguous_beds,
        [[],[]]
    )
    versions = versions.mix(HAPPY_SOMPY.out.versions.first())

    HAPPY_SOMPY.out.stats
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "sompy"], file) }
        .groupTuple()
        .set{ summary_reports }

    SPLIT_SOMPY_FEATURES(
        HAPPY_SOMPY.out.features
    )
    versions = versions.mix(SPLIT_SOMPY_FEATURES.out.versions)

    SPLIT_SOMPY_FEATURES.out.TP
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "sompy"], file) }
        .set{tp_vars_base}

    SPLIT_SOMPY_FEATURES.out.FP
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "sompy"], file) }
        .set{fp_vars}

    SPLIT_SOMPY_FEATURES.out.FN
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "sompy"], file) }
        .set{fn_vars}

    FILTER_VARIANTS(
        input_ch.map{meta, test, test_index, _truth, _truth_index, regions, target ->
            [ meta, test, test_index, regions, target ]}.join(SPLIT_SOMPY_FEATURES.out.TP)
        )
    versions = versions.mix(FILTER_VARIANTS.out.versions)

    // reheader tp_comp vcf files for tagged results
    FILTER_VARIANTS.out.csv
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "sompy"], file) }
        .set { tp_vars_comp }


    tagged_variants_csv = tagged_variants_csv
                            .mix(tp_vars_base)
                            .mix(fp_vars)
                            .mix(fn_vars)
                            .mix(tp_vars_comp)

    emit:
    summary_reports     // channel: [val(meta), reports]
    tagged_variants_csv // channel: [val(meta), csvs]
    versions            // channel: [versions.yml]
}
