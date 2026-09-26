test_that("collapseFactorTstats replaces matched factor columns with RMS summaries", {
    wide_df <- data.frame(
        phenotype_id = c("gene1", "gene2"),
        variant_id = c("snp1", "snp2"),
        t.factor_batch_B = c(3, 0),
        t.factor_batch_C = c(4, 12),
        t.populationYRI = c(5, 8),
        t.populationCEU = c(12, 15),
        t.genotype = c(2, -3),
        check.names = FALSE
    )

    collapsed <- collapseFactorTstats(wide_df)

    expect_false(any(c("t.factor_batch_B", "t.factor_batch_C",
                       "t.populationYRI", "t.populationCEU") %in%
                     names(collapsed)))
    expect_named(collapsed, c("phenotype_id", "variant_id", "t.genotype",
                              "t.batch", "t.population"))
    expect_equal(collapsed[["t.batch"]], c(sqrt((3^2 + 4^2) / 2),
                                           sqrt((0^2 + 12^2) / 2)))
    expect_equal(collapsed[["t.population"]], c(sqrt((5^2 + 12^2) / 2),
                                                sqrt((8^2 + 15^2) / 2)))
    expect_equal(collapsed[["t.genotype"]], wide_df[["t.genotype"]])
})

test_that("collapseFactorTstats leaves unmatched and singleton groups unchanged", {
    wide_df <- data.frame(
        phenotype_id = "gene1",
        variant_id = "snp1",
        t.factor_batch_B = 3,
        t.genotype = 2,
        check.names = FALSE
    )

    collapsed <- collapseFactorTstats(wide_df)

    expect_identical(collapsed, wide_df)
})
