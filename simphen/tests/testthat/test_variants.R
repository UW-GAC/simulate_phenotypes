context("variant tests")
library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(GENESIS)

test_that("variant effect", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.sel=1, verbose=FALSE)
    G <- as.vector(seqGetData(gds, "$dosage_alt"))

    eff1 <- variant_effect(G, h2=0.001, varComp=c(3,7))
    expect_equal(G > 0, eff1$Gbeta > 0)
    expect_equal(eff1$h2, 0.001)
    
    eff2 <- variant_effect(G, beta=eff1$beta, varComp=c(3,7))
    expect_equal(eff1$beta, eff2$beta)
    expect_equal(eff1$h2, eff2$h2)
    expect_equal(eff1$Gbeta, eff2$Gbeta)
    
    eff1 <- variant_effect(G, beta=0.5, varComp=c(3,7))
    expect_equal(G > 0, eff1$Gbeta > 0)
    expect_equal(eff1$beta, 0.5)
    
    eff2 <- variant_effect(G, h2=eff1$h2, varComp=c(3,7))
    expect_equal(eff1$beta, eff2$beta)
    expect_equal(eff1$h2, eff2$h2)
    expect_equal(eff1$Gbeta, eff2$Gbeta)
    
    seqClose(gds)
})


test_that("variant selection", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    
    af <- do.call(cbind, lapply(c("EAS_AF", "EUR_AF", "AFR_AF", "AMR_AF", "SAS_AF"), function(x) {
        tmp <- seqGetData(gds, paste0("annotation/info/", x))
        if ("data" %in% names(tmp)) tmp <- tmp$data
        tmp
    }))
    colnames(af) <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    # select variants in AFR but not EUR
    vs <- var_single_stratum(af[,"AFR"], af[,"EUR"])

    # select variants in both
    vm <- var_multiple_strata(af[,"AFR"], af[,"EUR"])

    expect_equal(length(intersect(vs, vm)), 0)
    
    seqClose(gds)
})


test_that("variant_genotypes", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    variant.id <- seqGetData(gds, "variant.id")
    sample.id <- seqGetData(gds, "sample.id")
    
    geno <- variant_genotypes(gds, variant.sel=1:10)
    expect_equal(colnames(geno), as.character(variant.id[1:10]))
    expect_equal(rownames(geno), sample.id)
    
    geno2 <- variant_genotypes(gds, variant.id=variant.id[1:10])
    expect_equal(geno, geno2)
    
    geno <- variant_genotypes(gds, variant.sel=1:10, sample.id=sample.id[1:10])
    expect_equal(rownames(geno), sample.id[1:10])
    
    seqClose(gds)
})


test_that("beta_pooled", {
    strata <- list(a=1:100, b=101:200)
    beta <- list(a=1, b=2)
    expect_equal(.beta_pooled(beta, strata), 1.5)
    
    beta <- list(a=c(1,1), b=c(2,3))
    expect_equal(.beta_pooled(beta, strata), c(1.5, 2))
})


test_that("variant assoc matches iterator method", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    data("1KG_phase3_annot")
    data("1KG_pcrelate_blockDiag_Matrix")
    varComp <- c(4,10)
    set.seed(7)
    outcome <- outcome_from_covMat_blocks(blockDiagMatrix1KG, varComp)
    dat <- annot
    dat$outcome <- outcome[dat$sample.id]
    
    var <- 1000
    beta <- 1.2
    geno <- variant_genotypes(gds, variant.sel=var)
    assoc <- variant_assoc(geno, beta=beta, varComp=varComp,
                           dat=dat, outcome="outcome", covars="Population",
                           cov.mat=blockDiagMatrix1KG)
    expect_equal(assoc$group, "all")

    eff <- variant_effect(G=geno, beta=beta, varComp=varComp)
    dat$outcome <- dat$outcome + eff$Gbeta
    nullmod <- fitNullModel(dat, outcome="outcome", covars="Population",
                            cov.mat=blockDiagMatrix1KG, 
                            verbose=FALSE)
    seqData <- SeqVarData(gds, sampleData=dat)
    seqSetFilter(gds, variant.sel=var, verbose=FALSE)
    iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)

    expect_equal(assoc$Est, assoc2$Est, tolerance=1e-6)
    
    seqClose(gds)
})


test_that("variant assoc", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    data("1KG_phase3_annot")
    data("1KG_pcrelate_blockDiag_Matrix")
    varComp <- c(4,10)
    set.seed(4)
    outcome <- outcome_from_covMat_blocks(blockDiagMatrix1KG, varComp)
    dat <- annot
    dat$outcome <- outcome[dat$sample.id]
    
    # select variants for testing
    set.seed(5)
    pruned <- SNPRelate::snpgdsLDpruning(gds, maf=0.05, verbose=FALSE)
    pruned.id <- unlist(pruned, use.names=FALSE)
    
    # select strata
    strata <- lapply(c("AFR", "EUR"), function(x) dat$sample.id[dat$Super.Population %in% x])
    names(strata) <- c("AFR", "EUR")
    
    beta <- list(AFR=1.2, EUR=1.5)
    geno <- variant_genotypes(gds, variant.id=pruned.id[1:3])
    assoc <- variant_assoc(geno[,1,drop=FALSE], strata=strata, beta=beta, varComp=varComp,
                           dat=dat, outcome="outcome", covars="Population",
                           cov.mat=blockDiagMatrix1KG)
    expect_true(all(abs(assoc$Est - assoc$beta) < assoc$Est.SE))
    expect_equal(assoc$group, c(names(strata), "pooled"))
    bp <- assoc$beta[assoc$group == "pooled"]
    expect_true(bp > 1.2)
    expect_true(bp < 1.5)
    
    # check fewer samples
    sm <- lapply(strata, function(x) x[1:200])
    assoc2 <- variant_assoc(geno[,1,drop=FALSE], strata=sm, beta=beta, varComp=varComp,
                            dat=dat, outcome="outcome", covars="Population",
                            cov.mat=blockDiagMatrix1KG)
    expect_true(all(assoc2$Est < assoc$Est))
    
    # rearrange sample.id
    set.seed(6)
    dat3 <- dat[sample(nrow(dat)),]
    assoc3 <- variant_assoc(geno[,1,drop=FALSE], strata=strata, beta=beta, varComp=varComp,
                            dat=dat3, outcome="outcome", covars="Population",
                            cov.mat=blockDiagMatrix1KG)
    expect_equal(assoc3$Est, assoc$Est)
    
    # multiple variants
    geno <- variant_genotypes(gds, variant.id=pruned.id[c(1,20,30)])
    assoc4 <- variant_assoc(geno, strata=strata, beta=beta, varComp=varComp,
                            dat=dat, outcome="outcome", covars="Population",
                            cov.mat=blockDiagMatrix1KG)
    expect_equal(nrow(assoc4), ncol(geno)*(length(strata)+1))
    hp <- assoc4$h2[assoc4$group == "pooled"]
    expect_true(all(!is.na(hp)))
    
    # use h2
    h2 <- list(AFR=assoc$h2[assoc$group == "AFR"], EUR=assoc$h2[assoc$group == "EUR"])
    assoc5 <- variant_assoc(geno[,1,drop=FALSE], strata=strata, h2=h2, varComp=varComp,
                           dat=dat, outcome="outcome", covars="Population",
                           cov.mat=blockDiagMatrix1KG)
    expect_equal(assoc5$Est, assoc$Est)
    
    seqClose(gds)
})
