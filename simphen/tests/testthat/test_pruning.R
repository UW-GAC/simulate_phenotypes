context("ld tests")
library(SeqArray)
library(SNPRelate)

test_that("ld pruning", {
    data("1KG_phase3_annot")
    super.pops <- unique(annot$Super.Population)
    strata <- lapply(super.pops, function(x) annot$sample.id[annot$Super.Population %in% x])
    names(strata) <- super.pops
    strata2 <- strata[c("AFR", "EUR")]
    
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    
    set.seed(7)
    pruned.all <- snpgdsLDpruning(gds, sample.id=unlist(strata2),
                                  maf=0.01, method="corr", ld.threshold=0.5,
                                  verbose=FALSE)
    pruned.all <- unlist(pruned.all, use.names=FALSE)
    
    set.seed(71)
    pruned.afr <- snpgdsLDpruning(gds, sample.id=strata2[["AFR"]], snp.id=pruned.all,
                                  remove.monosnp=FALSE, method="corr", ld.threshold=0.5,
                                  verbose=FALSE)
    pruned.afr <- unlist(pruned.afr, use.names=FALSE)
    
    set.seed(72)
    pruned.eur <- snpgdsLDpruning(gds, sample.id=strata2[["EUR"]], snp.id=pruned.afr,
                                  remove.monosnp=FALSE, method="corr", ld.threshold=0.5,
                                  verbose=FALSE)
    pruned <- unlist(pruned.eur, use.names=FALSE)
    
    chk <- ldprune_strata(gds, strata2, ld.threshold=0.5, seed=list(all=7, AFR=71, EUR=72))
    expect_equal(pruned, chk)
                         
    seqClose(gds)
})