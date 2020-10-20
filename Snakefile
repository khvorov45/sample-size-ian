rule install_deps:
    input:
        "renv.lock"
    output:
        touch(".deps-installed")
    shell:
        """Rscript -e 'renv::restore()'"""


rule sample_size_gmr:
    input:
        rules.install_deps.output,
        "gmr/gmr.R",
    output:
        "gmr/example-sample-diffs.pdf",
        "gmr/example-sample-titres.pdf",
        "gmr/results.csv",
        "gmr/summary.csv",
    shell:
        "Rscript gmr/gmr.R"

rule all:
    input:
        rules.install_deps.output,
        rules.sample_size_gmr.output,
