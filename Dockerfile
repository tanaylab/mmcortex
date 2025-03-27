# FROM bioconductor/bioconductor_docker:RELEASE_3_20
FROM r-base:4.4.2

# Install rpm dependencies
# RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libgtk2.0-dev libcairo2-dev libxt-dev xvfb xauth xfonts-base vim && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y  git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc zlib1g-dev libgtk2.0-dev libcairo2-dev libxt-dev xvfb && rm -rf /var/lib/apt/lists/*


RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("grid")'
RUN R -e 'install.packages("gridExtra")'
RUN R -e 'install.packages("pheatmap")'
RUN R -e 'install.packages("matrixStats")'
RUN R -e 'install.packages("umap")'
RUN R -e 'install.packages("Matrix")'
RUN R -e 'install.packages("princurve")'
RUN R -e 'install.packages("vioplot")'
RUN R -e 'install.packages("viridis")'
RUN R -e 'install.packages("cowplot")'
RUN R -e 'install.packages("patchwork")'
RUN R -e 'install.packages("ggrastr")'
RUN R -e 'install.packages("ggplotify")'
RUN R -e 'install.packages("Gviz")'
RUN R -e 'install.packages("gridExtra")'



RUN R -e 'if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}'
RUN R -e 'BiocManager::install(version = "3.20")'

RUN R -e 'remotes::install_github("tanaylab/tgstat")'
RUN R -e 'remotes::install_github("tanaylab/tglkmeans")'
RUN R -e 'BiocManager::install("tanaylab/metacell")'
RUN R -e 'BiocManager::install("ComplexHeatmap")'
RUN R -e 'remotes::install_github("tanaylab/mcatac")'
RUN R -e 'remotes::install_github("tanaylab/metacell.flow")'
RUN R -e 'remotes::install_github("tanaylab/misha")'
RUN R -e 'remotes::install_github("tanaylab/misha.ext")'
RUN R -e 'remotes::install_github("tanaylab/shaman")'
RUN R -e 'remotes::install_github("tanaylab/mcATAC")'
RUN R -e 'remotes::install_github("tanaylab/iceqream")'

RUN git clone https://yonatans2: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@github.com/tanaylab/mmcortex.git

WORKDIR /mmcortex

RUN R -e 'source("./scripts/download_data.r")'

RUN bash scripts/pl/pl_figures.sh

# Run R
# CMD ["R"]
