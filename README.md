This repo holds links and code to browse processed data and reproduce the figures in [Neural stem cell epigenomes and fate bias are temporally coordinated during corticogenesis](https://www.biorxiv.org/content/10.1101/2025.04.03.647013v1) by Shapira et al. (2025).

Figures can be reproduced in two ways:
1) Manually - downloading the repo and data separately, and installing required packages.
2) Using a docker container.

How to reproduce manually
1) Clone this repo
2) Install packages: run `Rscript scripts/install_requirements.r` 
3) In the repo's directory, download processed data using `wget  https://mmcortex.s3.eu-west-1.amazonaws.com/mmcortex.tar.gz`
4) Extract the data using `tar -xvzf mmcortex.tar.gz`.
5) Run `bash scripts/pl/pl_figures.sh`. The figures will be in `mmcortex/output/paper_figs`.

Using docker container
1) Install [docker](https://docs.docker.com/engine/install/).
2) Download [container](https://mmcortex.s3.eu-west-1.amazonaws.com/mmcortex_docker.tar).
3) Run container using the command `docker run -d --name <MY_DOCKER_NAME> -it mmcortex sh`. (choose container alias instead of <MY_DOCKER_NAME>)
4) Inside the container enter the command `bash scripts/pl/pl_figures.sh`.
5) Open a different terminal and change to the directory where you are running the container.
6) In the new terminal, copy output figures to your host environment using `docker cp <MY_DOCKER_NAME>:/mmcortex/output/paper_figs/ <MY_DESINATION_DIRECTORY>`
