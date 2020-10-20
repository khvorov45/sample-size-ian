rule install_deps:
    input:
        "renv.lock"
    output:
        touch(".deps-installed")
    shell:
        """Rscript -e 'renv::restore()'"""


rule gmr:
    input:
        rules.install_deps.output,
        "gmr/gmr.R",
    output:
        "gmr/example-sample-diffs.pdf",
        "gmr/example-sample-titres.pdf",
        "gmr/results.csv",
        "gmr/summary-gmr.csv",
    shell:
        "Rscript gmr/gmr.R"

rule seroconv:
    input:
        rules.install_deps.output,
        "seroconv/seroconv.R",
    output:
        "seroconv/example-table.csv",
        "seroconv/results.csv",
        "seroconv/summary-seroconv.csv",
    shell:
        "Rscript seroconv/seroconv.R"

rule all:
    input:
        rules.install_deps.output,
        rules.gmr.output,
        rules.seroconv.output,
