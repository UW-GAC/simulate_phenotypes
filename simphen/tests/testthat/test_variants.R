context("variant tests")
library(SeqArray)

test_that("variant effect", {
    gdsfmt::showfile.gds(closeall=TRUE)
    gdsfile <- system.file("extdata", "1KG_phase3_chr1_SNVsubset.gds", package="simphen")
    gds <- seqOpen(gdsfile)
    seqSetFilter(gds, variant.sel=1, verbose=FALSE)
    G <- as.vector(seqGetData(gds, "$dosage_alt"))

    outcome <- rnorm(length(G))
    outc_eff <- add_variant_effect(outcome, G, h2=0.001)
    expect_true(all(outc_eff - outcome >= 0))
    
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
