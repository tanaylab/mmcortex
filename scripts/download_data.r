# script for downloading data for reproducing figures of the mmcortex paper

options(timeout = 100000)

download.file("https://mmcortex.s3.eu-west-1.amazonaws.com/mmcortex.tar.gz", "mmcortex.tar.gz")

system("tar -xvzf mmcortex.tar.gz")

folder_names <- c('figs', 
                  'output/paper_figs')

for (folder in folder_names) {
  if(!dir.exists(folder)) {
    dir.create(folder)
  }
  print(list.files(folder))
}
