# script for downloading data for reproducing figures of the mmcortex paper

options(timeout = 100000)

download.file("https://mmcortex.s3.eu-west-1.amazonaws.com/mmcortex.tar.gz", "mmcortex.tar.gz")

system("tar -xvzf mmcortex.tar.gz")

# file.remove("mmcortex.tar.gz")

folder_names <- c('figs', 
#                   'scripts', 
#                   'output', 
#                   'scdb',
#                   'output/metacell_model', 
#                   'output/metacell_flow', 
#                   'output/mcatac',
#                   'output/methylation',
#                   'output/hic',
#                   'output/sequence_modeling',
#                   'output/MPRA',
                  'output/paper_figs')

for (folder in folder_names) {
  if(!dir.exists(folder)) {
    dir.create(folder)
  }
  print(list.files(folder))
}
