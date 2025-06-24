FROM r-base:latest

RUN apt-get update && \
    apt-get install -y libssl-dev libcurl4-openssl-dev git cmake libprotobuf-dev protobuf-compiler

WORKDIR /usr/local/bin

RUN wget -O hmmcopy_utils.zip https://codeload.github.com/shahcompbio/hmmcopy_utils/zip/refs/heads/master && \
    unzip hmmcopy_utils.zip && \
    git clone https://github.com/broadinstitute/ichorCNA.git && \
    ln -s /usr/local/bin/ichorCNA/scripts/runIchorCNA.R /usr/local/bin/runIchorCNA.R  

RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('GenomicRanges', 'HMMcopy', 'GenomeInfoDb', 'plyr', 'optparse'), update = TRUE, ask = FALSE)"

WORKDIR /usr/local/bin/hmmcopy_utils-master
RUN cmake . && \ 
    make && \
    ln -s /usr/local/bin/hmmcopy_utils-master/bin/readCounter /usr/local/bin/readCounter

RUN R CMD INSTALL /usr/local/bin/ichorCNA  

CMD ["/init"]
