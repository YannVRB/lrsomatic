process DMR_CHECK_MODIFICATIONS {
    tag "${meta.id}"
    label 'process_single'

    // Same image as DMR_HAPLOTYPE_REGIONS (it only needs zcat and awk)
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f1/f1a9e30012e1b41baf9acd1ff94e01161138d8aa17f4e97aa32f2dc4effafcd1/data'
        : 'community.wave.seqera.io/library/pybedtools_bedtools_htslib_pip_pypints:39699b96998ec5f6'}"

    input:
    tuple val(meta), path(bedmethyl)

    output:
    tuple val(meta), path(bedmethyl), emit: checked

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'check_modifications.sh'

    stub:
    """
    true
    """
}
