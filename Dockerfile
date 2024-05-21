FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN Rscript -e 'BiocManager::install("GENESIS")'
RUN Rscript -e 'remotes::install_cran(c("argparser", "doParallel"))'
RUN Rscript -e 'remotes::install_github("UW-GAC/simulate_phenotypes/simphen")'
RUN cd /usr/local  && \
    git clone https://github.com/UW-GAC/simulate_phenotypes.git
