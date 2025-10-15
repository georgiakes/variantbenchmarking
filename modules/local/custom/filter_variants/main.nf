process FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/04/04b086608e967caf920b0af1b20c388933a2e600bb657fc4c4a119aeb15e6b6a/data':
        'community.wave.seqera.io/library/pysam_pandas:baaaf9c9ffb0e943' }"

    input:
    tuple val(meta), path(vcf), path(index), path(regions), path(targets), path(csv)

    output:
    tuple val(meta), path("*intersect.csv"), emit: csv
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = regions ? "--bed ${regions}" : ''
    def target = targets ? "--targets ${targets}" : ''

    """
    filter_variants.py \\
        --csv $csv \\
        --vcf $vcf \\
        $bed \\
        $target \\
        --out ${prefix}.intersect.csv
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.intersect.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
