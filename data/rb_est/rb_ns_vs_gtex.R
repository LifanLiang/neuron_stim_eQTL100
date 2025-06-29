library(data.table)

setwd("/project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/final_eqtl_08222023/")

eqtlfs <- list.files(pattern="trans2.5e5_eqtl.txt$")
eqtls <- list()
for(f in eqtlfs) {
  n <- substr(f, 1, nchar(f)-37)
  #eqtls[[n]] <- readRDS(f)$cis$eqtls
  eqtls[[n]] <- fread(f)
}

ref_pairs <- readRDS("../rb_est/ref_pairs.rds")

gtex.store <- "/scratch/midway3/lifanl/"
gtexfs <- list.files(gtex.store, pattern="Brain_")

snpinfo <- fread("../../snpinfo_100lines.txt")
snpinfo$variant_id <- paste(snpinfo$chr, snpinfo$pos, snpinfo$ref, snpinfo$alt, "b38", sep="_")

#cerebellum <- fread("../../GTEx_Analysis_v8_eQTL/Brain_Cerebellum.v8.signif_variant_gene_pairs.txt.gz")
#geneinfo <- fread("../../GTEx_Analysis_v8_eQTL/Brain_Cerebellum.v8.egenes.txt.gz")
#cere1 <- merge(cerebellum, geneinfo[,1:2], by="gene_id", all.x=T)
#cere2 <- merge(cere1, snpinfo[,c("rsid","variant_id")], by="variant_id", all=F)

## Estimate rb between our eQTL and GTEx Brain tissues
rb_est <- function(eqtl_ns, eqtl_gtex, ref_pairs) {
    null_ns <- eqtl_ns[gene %in% ref_pairs$gene & `p-value`>0.01,1:3]
    null_gtex <- eqtl_gtex[gene_id %in% ref_pairs$gene_id & `pval_nominal`>0.01,
                           c("slope","gene_name","rsid")]
    null.sumstats <- merge(null_ns, null_gtex, by.x=c("SNP","gene"), 
                           by.y=c("rsid","gene_name"), all=F)
    dtCor <- null.sumstats[, .(mCor = cor(beta,slope)), by=gene]
    re <- mean(dtCor$mCor, na.rm=T)
    
    pair1 <- merge(ref_pairs[,1:2], eqtl_ns, by=c("gene","SNP"), all=F)
    pair1$se <- pair1$beta / pair1[,"t-stat"]
    pair2 <- merge(ref_pairs[,3:4], eqtl_gtex, by=c("gene_id","variant_id"), all=F)
    #pair2$se <- pair2$beta / pair2[,"t-stat"]

    var.e1 <- mean(pair1$se^2)
    var.e2 <- mean(pair2$slope_se^2)

    var.b1 <- var(pair1$beta)
    var.b2 <- var(pair2$slope)

    pairs <- merge(pair1, pair2, by.x=c("gene","SNP"), 
                   by.y=c("gene_name", "rsid"), all=F)
    covar.bhat <- cov(pairs$beta, pairs$slope)

    rb <- (covar.bhat - re*sqrt(var.e1*var.e2))/sqrt(var.b1-var.e1)/sqrt(var.b2-var.e2)
    rb
}

## Standard error for rb between our eqtl and GTEx
rb_se <- function(eqtl_ns, eqtl_gtex, ref_pairs) {
    null_ns <- eqtl_ns[gene %in% ref_pairs$gene & `p-value`>0.01,1:3]
    null_gtex <- eqtl_gtex[gene_id %in% ref_pairs$gene_id & `pval_nominal`>0.01,
                           c("slope","gene_name","rsid")]
    null.sumstats <- merge(null_ns, null_gtex, by.x=c("SNP","gene"), 
                           by.y=c("rsid","gene_name"), all=F)
    dtCor <- null.sumstats[, .(mCor = cor(beta,slope)), by=gene]
    
    pair1 <- merge(ref_pairs[,1:2], eqtl_ns, by=c("gene","SNP"), all=F)
    pair1$se <- pair1$beta / pair1[,"t-stat"]
    pair2 <- merge(ref_pairs[,3:4], eqtl_gtex, by=c("gene_id","variant_id"), all=F)
    
    rbs <- list()
    for(g in ref_pairs$gene) { #Jack-knife approximation
        re <- mean(dtCor$mCor[dtCor$gene!=g], na.rm=T)
    
        var.e1 <- mean(pair1$se[pair1$gene!=g]^2)
        var.e2 <- mean(pair2$slope_se[pair2$gene_name!=g]^2)

        var.b1 <- var(pair1$beta[pair1$gene!=g])
        var.b2 <- var(pair2$slope[pair2$gene_name!=g])

        pairs <- merge(pair1[pair1$gene!=g], pair2[pair2$gene_name!=g], 
                       by.x=c("gene","SNP"), by.y=c("gene_name", "rsid"), all=F)
        covar.bhat <- cov(pairs$beta, pairs$slope)
        #covar.bhat <- cov(pair1$beta[pair1$gene!=g], pair2$beta[pair2$gene!=g])

        rbs[[g]] <- (covar.bhat - re*sqrt(var.e1*var.e2))/sqrt(var.b1-var.e1)/sqrt(var.b2-var.e2)
    }
    rb1 <- unlist(rbs)
    sqrt((length(rb1)-1)/length(rb1)*sum((rb1-mean(rb1))^2)) # Standard error of rb
}


res <- matrix(nrow=length(eqtls),ncol=length(gtexfs))
rownames(res) <- names(eqtls)
colnames(res) <- substr(gtexfs,7,nchar(gtexfs)-16)
res.se <- res
for(f in gtexfs) {
    n <- substr(f,7,nchar(f)-16)
    temp <- fread(paste0(gtex.store,f))
    geneinfo <- fread(paste0("../../GTEx_Analysis_v8_eQTL/Brain_",n,".v8.egenes.txt.gz"))
    temp1 <- merge(temp, geneinfo[,1:2], by="gene_id", all.x=T)
    temp2 <- merge(temp1, snpinfo[,c("rsid","variant_id")], by="variant_id", all=F)
    for(e in names(eqtls)) {
        cat("Estimate rb and its standard error between ",e," and ", n)
        res[e, n] <- rb_est(eqtls[[e]], temp2, ref_pairs)
        res.se[e, n] <- rb_se(eqtls[[e]], temp2, ref_pairs)
    }
}
saveRDS(list(rb=res,rb.se=res.se),"/project/xinhe/lifan/neuron_stim/mateqtl_100lines_output/rb_est/rb_ns_vs_gtex.rds")
