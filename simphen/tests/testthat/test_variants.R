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
        seqGetData(gds, paste0("annotation/info/", x))$data
    }))
    colnames(af) <- c("EAS", "EUR", "AFR", "AMR", "SAS")

    # select variants in AFR but not EUR
    vs <- var_single_stratum(af[,"AFR"], af[,"EUR"])

    # select variants in both
    vm <- var_multiple_strata(af[,"AFR"], af[,"EUR"])

    expect_equal(length(intersect(vs, vm)), 0)
    
    seqClose(gds)
})


test_that("variant assoc", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    data("1KG_pcrelate_blockDiag_Matrix")
    varComp <- c(4,10)
    set.seed(4)
    outcome <- outcome_from_covMat_blocks(blockDiagMatrix1KG, varComp)
    set.seed(5)
    covars <- data.frame(x=rnorm(ncol(blockDiagMatrix1KG)))
    dat <- Biobase::AnnotatedDataFrame(cbind(data.frame(outcome), covars, sample.id=names(outcome), stringsAsFactors=FALSE))

    beta <- 1.2
    assoc <- variant_assoc(variant.sel=18, beta=beta, varComp=varComp,
                           gdsobj=gds, dat=dat,
                           outcome="outcome", covars="x",
                           cov.mat=blockDiagMatrix1KG)
    expect_equal(assoc$Est, beta, tolerance=0.01)

    # check fewer samples
    assoc2 <- variant_assoc(variant.sel=18, beta=beta, varComp=varComp,
                           gdsobj=gds, dat=dat[1:1000,],
                           outcome="outcome", covars="x",
                           cov.mat=blockDiagMatrix1KG)
    expect_true(assoc2$Est < assoc$Est)

    # rearrange sample.id
    set.seed(6)
    dat3 <- dat[sample(nrow(dat)),]
    assoc3 <- variant_assoc(variant.sel=18, beta=beta, varComp=varComp,
                            gdsobj=gds, dat=dat3,
                            outcome="outcome", covars="x",
                            cov.mat=blockDiagMatrix1KG)
    expect_equal(assoc3$Est, assoc$Est)
    
    seqClose(gds)
})



test_that("variant assoc matches iterator method", {
    gdsfmt::showfile.gds(closeall=TRUE, verbose=FALSE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    data("1KG_pcrelate_blockDiag_Matrix")
    varComp <- c(4,10)
    set.seed(7)
    outcome <- outcome_from_covMat_blocks(blockDiagMatrix1KG, varComp)
    set.seed(8)
    covars <- data.frame(x=rnorm(ncol(blockDiagMatrix1KG)))
    dat <- Biobase::AnnotatedDataFrame(cbind(data.frame(outcome), covars, sample.id=names(outcome), stringsAsFactors=FALSE))

    var <- 132
    beta <- 1.2
    assoc <- variant_assoc(variant.sel=var, beta=beta, varComp=varComp,
                           gdsobj=gds, dat=dat,
                           outcome="outcome", covars="x",
                           cov.mat=blockDiagMatrix1KG)

    seqSetFilter(gds, variant.sel=var, verbose=FALSE)
    geno <- as.vector(seqGetData(gds, "$dosage_alt"))
    eff <- variant_effect(G=geno, beta=beta, varComp=varComp)
    dat <- dat[seqGetData(gds, "sample.id"),]
    dat$outcome <- dat$outcome + eff$Gbeta
    nullmod <- fitNullModel(dat, outcome="outcome", covars="x",
                            cov.mat=blockDiagMatrix1KG, 
                            verbose=FALSE)
    seqData <- SeqVarData(gds, sampleData=dat)
    seqSetFilter(gds, variant.sel=var, verbose=FALSE)
    iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
    assoc2 <- assocTestSingle(iterator, nullmod, verbose=FALSE)

    expect_equal(assoc$Est, assoc2$Est)
    
    seqClose(gds)
})
