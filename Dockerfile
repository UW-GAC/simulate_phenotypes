FROM uwgac/topmed-master:2.12.0

USER root
RUN cd /usr/local  && \
    git clone https://github.com/UW-GAC/simulate_phenotypes.git && \
    cd simulate_phenotypes && \
    R CMD INSTALL simphen
