FROM rocker/tidyverse:latest

RUN apt-get update && \
    apt-get -y install libbz2-dev liblzma-dev libproj-dev libgdal-dev cmake libcurl4-openssl-dev libssl-dev libxml2-dev

RUN R -e 'install.packages("librarian")' && \
    R -e 'BiocManager::install("Biobase")'
RUN R -e 'librarian::shelf(caret, CGHcall, data.table, DBI, DescTools, devtools, doMC, dplyr, EnsDb.Hsapiens.v75, \
    flexmix, GenomicRanges, ggalt, ggplot2, ggpubr, ggrepel, gridExtra, hrbrthemes, ks, magrittr, NMF, pheatmap, plyr, \
    purrr, readr, RMySQL, stringr, tidyr, viridis, YAPSA, glue, doMC, EnsDb.Hsapiens.v86, EnsDb.Mmusculus.v79)'
RUN R -e 'librarian::shelf(QDNAseq)'

WORKDIR /usr/local/lib/R/site-library/
RUN R -e "BiocManager::install('YAPSA', update = TRUE, ask = FALSE)"
RUN R -e "devtools::install_github('stephenturner/annotables', update = TRUE, ask = FALSE, dependencies = TRUE)" 
RUN R -e "devtools::install_github('Huntsmanlab/utanos', update = TRUE, ask = FALSE, dependencies = TRUE)"

CMD ["/init"]
