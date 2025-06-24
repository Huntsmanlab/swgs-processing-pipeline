FROM rocker/tidyverse:latest

RUN apt-get update

RUN R -e 'install.packages("AzureStor")'
RUN R -e 'install.packages("processx")'
ENV PATH="$PATH"
WORKDIR /usr/bin
RUN wget https://aka.ms/downloadazcopy-v10-linux && \
    tar -xvzf downloadazcopy-v10-linux && \
    ln -s azcopy_linux_amd64_10.28.0/azcopy azcopy && \
    rm -rf azcopy_linux_amd64_10.26.0 downloadazcopy-v10-linux

CMD ["/init"]
