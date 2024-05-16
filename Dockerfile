FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN Rscript -e 'BiocManager::install("GENESIS")'
RUN Rscript -e 'remotes::install_cran("doParallel")'
RUN Rscript -e 'remotes::install_github("UW-GAC/simulate_phenotypes/simphen")'
