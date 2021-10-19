FROM uwgac/topmed-master:2.8.0

USER root
RUN cd /usr/local  && \
    git clone -b v0.2.2 https://github.com/UW-GAC/simulate_phenotypes.git && \
    cd simulate_phenotypes && \
    R CMD INSTALL simphen && \
    Rscript -e 'install.packages("doParallel", repos="https://cloud.r-project.org")'
