library(pheatmap)

setwd("/project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/rb_est")

rb.gtex <- readRDS("rb_ns_vs_gtex.rds")
cell_type.idx <- c(1,4,7,2,5,8,3,6,9)
dat <- rb.gtex$rb[cell_type.idx,]
#dat <- rb.gtex$rb
dat <- dat[,order(colMeans(dat),decreasing=T)]
dat <- dat[,colnames(dat)!="Cortex"] # Cortex top eQTL was used to construct the reference set
pheatmap(dat, cluster_rows = F, cluster_cols = F, angle_col="315")

ns_rb <- read.table("neuron_stim_rb_ref_cortex.txt")
dat1 <- ns_rb[cell_type.idx,cell_type.idx]
diag(dat1) <- NA
pheatmap(dat1, cluster_rows = F, cluster_cols = F, angle_col="315",
         color = colorRampPalette(c("white", "red"))(100))
