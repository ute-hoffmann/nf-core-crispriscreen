process PREPARE_LIBRARY {
    tag "$prepare_libary"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::r-base=4.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' :
        'quay.io/biocontainers/mulled-v2-8849acf39a43cdd6c839a369a74c0adc823e2f91:ab110436faf952a33575c64dd74615a84011450b-0' }"

    input:
    path fasta

    output:
    path '*.saf'       , emit: annotation
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/crispriscreen/bin/
    """
    prepare_library.R \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
