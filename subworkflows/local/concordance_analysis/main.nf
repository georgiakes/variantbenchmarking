//
// CONCORDANCE_ANALYSIS: SUBWORKFLOW FOR CONCORDANCE ANALYSIS BETWEEN BASE AND TEST VCFS
//

include { GATK4_CONCORDANCE                 } from '../../../modules/nf-core/gatk4/concordance'
include { PICARD_CREATESEQUENCEDICTIONARY   } from '../../../modules/nf-core/picard/createsequencedictionary'

workflow CONCORDANCE_ANALYSIS {
    take:
    input_ch
    bed_ch
    fasta_ch
    fai_ch
    dictionary

    main:

    versions        = Channel.empty()

    //prepare dict file for liftover of vcf files
    if (!params.dictionary){
        PICARD_CREATESEQUENCEDICTIONARY(
            fasta_ch
        )
        dictionary = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        versions = versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
    }

    ch_pairs = input_ch
        .map { meta, vcf1, tbi1 -> 
            // Simplify to just keep the id and files
            [meta.id, vcf1, tbi1] 
        }
        .toList()
        .flatMap { items ->
            def result = []
            
            // Generate pairwise combinations
            for (int i = 0; i < items.size(); i++) {
                for (int j = i + 1; j < items.size(); j++) {
                    def left = items[i]   // [id1, vcf1, tbi1]
                    def right = items[j]  // [id2, vcf2, tbi2]
                    
                    // Create new metadata with combined IDs
                    def combinedMeta = [id: "${left[0]}-${right[0]}"]
                    
                    result << [
                        combinedMeta,     // [id: "test7-test6"]
                        left[1], left[2], // vcf1, tbi1 from first sample
                        right[1], right[2] // vcf2, tbi2 from second sample
                    ]
                }
            }
            return result
        }
            
    GATK4_CONCORDANCE(
        ch_pairs,
        bed_ch,
        fasta_ch,
        fai_ch,
        dictionary
    )
    versions = versions.mix(GATK4_CONCORDANCE.out.versions)

    emit:
    versions     // channel: [val(meta), versions.yml]
}
