FROM r-base:latest

RUN apt-get update && \
    apt-get install -y libssl-dev libcurl4-openssl-dev

RUN R -e "install.packages('BiocManager')" 
RUN R -e "BiocManager::install(c('QDNAseq', 'ACE', 'BSgenome.Mmusculus.UCSC.mm10', 'BSgenome.Hsapiens.UCSC.hg19', 'BSgenome.Hsapiens.UCSC.hg38', 'dplyr', 'Biobase', 'ggplot2', 'stringr', 'future', 'CGHcall', 'glue'), update = TRUE, ask = FALSE)"

CMD ["/init"]
