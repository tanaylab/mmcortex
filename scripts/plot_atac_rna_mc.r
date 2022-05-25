gene = 'Hes6'
st = 'Cortex or Hippocampal Neuroblast'

library(metacell)
library(tidyverse)
wd = '/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex'
scdb_init(file.path(wd, 'scdb'), force_reinit = T)
scfigs_init(file.path(wd, 'figs'))
mc_prom = scdb_mc('prom')
mc_all = scdb_mc('all')
mc_md = vroom::vroom('/net/mraid14/export/tgdata/users/yonshap/proj/mmcortex/BonevCollab/mc_metadata.tsv')
df3 = left_join(mc_md %>% group_by(inferred_subtype) %>% select(color) %>% unique, 
                mc_md %>% group_by(inferred_subtype) %>% count,
                by = 'inferred_subtype') %>% arrange(desc(n))
mc_md$inferred_subtype = factor(mc_md$inferred_subtype, levels = df3$inferred_subtype)
mc_md$color = factor(mc_md$color, levels = df3$color)


st_inds = which(mc_md$inferred_subtype == st)
new_df = data.frame(cbind(mc_md$mean_day[st_inds], 
                          mc_prom@mc_fp[grep(gene, rownames(mc_prom@mc_fp)),st_inds],
                         mc_all@mc_fp[grep(gene, rownames(mc_all@mc_fp)),st_inds]))
colnames(new_df) = c('day', 'atac', 'rna')
new_df[,-1] = apply(new_df[,-1], 2, function(x) x/max(x))


df_pvt = pivot_longer(mutate(new_df, atac = atac/max(atac), rna = rna/max(rna)), cols = c('atac', 'rna')) %>% arrange(name, day)

p = ggplot(df_pvt, aes(x=day, y = value, color = name)) + geom_point() + geom_smooth() + ggtitle(paste(gene, '-', st))

ggsave(file.path(wd, 'figs', paste0(gene, '_', st, '_mc_atac_rna.png')), plot = p, width = 6, height = 6)