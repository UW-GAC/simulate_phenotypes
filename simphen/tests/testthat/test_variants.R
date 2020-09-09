context("variant tests")
library(SeqArray)

test_that("variant effect", {
    gdsfmt::showfile.gds(closeall=TRUE)
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
    gdsfmt::showfile.gds(closeall=TRUE)
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
