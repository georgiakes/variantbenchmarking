process GATK4_CONCORDANCE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    path(bed)
    tuple val(meta4), path(fasta)
    tuple val(meta5), path(fai)
    tuple val(meta6), path(dict)

    output:
    tuple val(meta), path('*.tsv')     ,    emit: summary
    tuple val(meta), path("*.tpfn.vcf"),    emit: tpfn
    tuple val(meta), path("*.tpfp.vcf"),    emit: tpfp
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = bed ? "--intervals $bed" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK MergeVcfs] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        Concordance \\
        -R $fasta \\
        -eval $vcf \\
        --truth $truth_vcf \\
        --summary ${prefix}.summary.tsv \\
        -tpfn ${prefix}.tpfn.vcf \\
        -tpfp ${prefix}.tpfp.vcf \\
        $intervals \\
        --tmp-dir . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.summary.tsv 
    touch ${prefix}.tpfn.vcf
    touch ${prefix}.tpfp.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
