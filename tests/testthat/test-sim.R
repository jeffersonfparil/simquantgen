test_that(
    "SIMULATIONS", {
        print("SIMULATIONS:")
        l = 500
        ploidy=4
        n_alleles=3
        n_effects = 10
        h2=0.75
        pheno_reps=3
        env_factor_levels=c(1, 2, 3, 2, 5)
        G = fn_simulate_genotypes(l=l, ploidy=ploidy, n_alleles=n_alleles)
        y = fn_simulate_phenotypes(G=G, n_alleles=n_alleles, n_effects=n_effects, h2=h2)$Y
        y_complex = fn_simulate_phenotypes(G=G, n_alleles=n_alleles, n_effects=n_effects, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, h2=h2)$Y
        Y = fn_simulate_phenotypes(G=G, n_alleles=n_alleles, n_effects=n_effects, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, h2=h2, pheno_reps=pheno_reps)$Y
        allele_freq_sum_per_locus = rowSums(matrix(colMeans(G), ncol=(n_alleles-1), byrow=TRUE)); # which(allele_freq_sum_per_locus > 1)
        list_df_CORR = fn_simulate_gxe(G=G, n_alleles=n_alleles, env_factor_levels=env_factor_levels, env_factor_effects_sd=0.75, n_reps=2)
        list_y_complete_y_X_Z_D = fn_simulate_sparse_mixed_model_data(G=G, df_gxe=list_df_CORR$df, frac_gen_missing=0.25)
        expect_equal(cor(G[,1], G[,n_alleles]) > cor(G[,1], G[,100]), TRUE)
        expect_equal(sum(allele_freq_sum_per_locus <= 1), l)
        expect_equal(length(y), nrow(G))
        expect_equal(length(y_complex), nrow(G))
        expect_equal(nrow(Y), nrow(G))
        expect_equal(ncol(Y), pheno_reps)
        expect_equal(ncol(list_df_CORR$df), 4+length(env_factor_levels))
        expect_equal(nrow(list_df_CORR$df), 100*prod(env_factor_levels)*2)
        expect_equal(ncol(list_df_CORR$CORR), prod(env_factor_levels))
        expect_equal(length(list_y_complete_y_X_Z_D$y_complete), length(list_y_complete_y_X_Z_D$y))
        expect_equal(nrow(list_df_CORR$df), length(list_y_complete_y_X_Z_D$y))
        expect_equal(sum(is.na(list_y_complete_y_X_Z_D$y)), floor(length(list_y_complete_y_X_Z_D$y)*0.25))
        expect_equal(nrow(list_y_complete_y_X_Z_D$y), nrow(list_y_complete_y_X_Z_D$X))
        expect_equal(nrow(list_y_complete_y_X_Z_D$y), nrow(list_y_complete_y_X_Z_D$Z))
        expect_equal(ncol(list_y_complete_y_X_Z_D$Z), nrow(list_y_complete_y_X_Z_D$D))
        expect_equal(nrow(list_y_complete_y_X_Z_D$D), ncol(list_y_complete_y_X_Z_D$D))
    }
)
