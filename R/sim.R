#' Simulate genotype data with some linkage disequilibrium between neighbouring loci
#' 
#' @param n: number of samples (individuals or pools)
#' @param l: number of loci
#' @param ploidy: ploidy level of individual samples or the number of individuals multiplied by their ploidy to simulate pools
#' @param n_alleles: number of alleles per locus
#' @param min_allele_freq: minimum minor allele frequency
#' @param n_chr: number of chromosomes
#' @param max_pos: total length of the genome
#' @param dist_bp_at_50perc_r2: distance in bases at which the estimated linkage between loci at both ends is at 50%
#' @param n_threads: number of computing cores or threads to use in parallel simulation of genotypes
#' @param verbose: show simulation messages?
#' @param show_correlation_heatmap: show correlation heatmap?
#' @returns genotype matrix with $n$ rows, $l x (n_alleles-1)$ columns, and named rows and columns
#' @examples
#' G_1 = fn_simulate_genotypes()
#' G_2 = fn_simulate_genotypes(n=150, l=10000, ploidy=4, n_alleles=3, verbose=TRUE)
#' G_3 = fn_simulate_genotypes(n=200, l=2e5, ploidy=4, n_alleles=3, n_chr=10, max_pos=2.2e9, dist_bp_at_50perc_r2=5e6, n_threads=32, verbose=TRUE)
#' @export
fn_simulate_genotypes = function(n=100, l=500, ploidy=2, n_alleles=2, min_allele_freq=0.01, n_chr=5, max_pos=135e6, dist_bp_at_50perc_r2=5e6, n_threads=2, verbose=FALSE, show_correlation_heatmap=FALSE) {
    # n=100
    # l=500
    # ploidy=2
    # n_alleles=2
    # min_allele_freq=0.01
    # n_chr=5
    # max_pos=135e6
    # dist_bp_at_50perc_r2=5e6
    # n_threads=2
    # verbose=FALSE
    # show_correlation_heatmap = TRUE
    ### Define loci positions
    vec_chr_lengths = rep(floor(l/n_chr), n_chr)
    vec_chr = rep(paste0("chr_", 1:n_chr), each=floor(l/n_chr))
    if (sum(vec_chr_lengths) < l) {
        n_missing = l - sum(vec_chr_lengths)
        vec_chr_lengths[n_chr] = vec_chr_lengths[n_chr] + n_missing
        vec_chr = c(vec_chr, rep(tail(vec_chr, n=1), each=n_missing))
    }
    vec_pos = c()
    for (i in 1:n_chr) {
        vec_pos = c(vec_pos, sort(sample.int(n=floor(max_pos/n_chr), size=vec_chr_lengths[i], replace=FALSE)))
    }
    ### Define the rate at which linkage decreases as a function of distance using a generalised logistic function
    if (verbose==TRUE) {
        print("Modelling linkage decay as a function of distance using a generalised logistic function.")
    }
    x = seq(0, 1, length=100)
    sigma = 1 / (dist_bp_at_50perc_r2 / floor(max_pos/n_chr))
    if (verbose==TRUE) {
        txtplot::txtplot(x, 1 / exp(sigma*x), xlab="Distance", ylab="Covariance")
    }
    ### Define the variance covariance matrix between loci
    if (verbose==TRUE) {
        print("Sampling allele frequencies from a multivariate distribution with variance-covariance matrix simulation LD within chromosome.")
    }
    m = vec_chr_lengths[1]
    x = seq(from=0, to=1, length=m)
    S = matrix(0, nrow=m, ncol=m)
    for (i in 1:m) {
        r = 0.25 / exp(sigma*x[1:(m-i+1)])
        S[i, i:m] = r
        S[i:m, i] = r
    }
    if (verbose==TRUE) {
        pb = txtProgressBar(min=0, max=n_chr, style=3)
    }
    for (i in 1:(n_chr-1)) {
        if (i==1) {
            Z = MASS::mvrnorm(n=n, mu=runif(n=m, min=min_allele_freq, max=(1-min_allele_freq)), Sigma=S)
        } else {
            Z = cbind(Z, MASS::mvrnorm(n=n, mu=runif(n=m, min=min_allele_freq, max=(1-min_allele_freq)), Sigma=S))
        }
        if (verbose==TRUE) {
            setTxtProgressBar(pb, i)
        }
    }
    ### Last chromosome
    m = vec_chr_lengths[n_chr]
    if (vec_chr_lengths[1] != m) {
        x = seq(from=0, to=1, length=m)
        S = matrix(0, nrow=m, ncol=m)
        for (i in 1:m) {
            r = 1 / exp(sigma*x[1:(m-i+1)])
            S[i, i:m] = r
            S[i:m, i] = r
        }
    }
    Z = cbind(Z, MASS::mvrnorm(n=n, mu=runif(n=m, min=min_allele_freq, max=(1-min_allele_freq)), Sigma=S))
    if (verbose==TRUE) {
        setTxtProgressBar(pb, n_chr)
        close(pb)
    }
    ### Scale the genotype data between maf and 1-maf
    if (verbose==TRUE) {
        print("Calling genotypes accounting for ploidy (multi-threaded across $n$ entries)...")
    }
    Z[Z < 0.0] = 0.0
    Z[Z > 1.0] = 1.0
    ### Define the expected genotype levels according to the ploidy
    ### Setting equal probabilities along the 0 to 1 range for each genotype level using quantiles
    geno_levels = c(0:ploidy)/ploidy
    n_genotype_levels = ploidy + 1
    mean_Z = mean(Z)
    sd_Z = sd(Z)
    q_steps = qnorm(p=c(1:(n_genotype_levels-1))/n_genotype_levels, mean=mean_Z, sd=sd_Z, lower.tail=TRUE)
    q_steps = c(0, q_steps)
    ### Sample alleles counts and convert them into allele frequencies
    clusters = parallel::makeCluster(n_threads, type="FORK")
	doParallel::registerDoParallel(clusters)
    `%dopar%` = foreach::`%dopar%`
    G = foreach::foreach(i=1:n, .combine=cbind, .inorder=TRUE) %dopar% {
		sapply(i, function(i) {
            g = rep(0, each=(l*(n_alleles-1)))
            for (j in 1:l) {
                for (k in 0:(n_alleles-2)) {
                    # i=1; j=1; k=1
                    a = (j-1)*(n_alleles-1) + 1
                    b = a + k
                    ### Skip the rest of the alleles if we have already selected the two alleles
                    if (sum(g[a:b]) == 1) {
                        break
                    }
                    if (k == 0) {
                        ### Set the first allele
                        idx = tail(which(Z[i, j] >= q_steps), n=1)
                        g[b] = geno_levels[idx]
                    } else {
                        ### Find the other allele where we expect a maximum of 2 non-zero alleles since we are still dealing with 2 sets of chromosomes split across gametes regardless of the even ploidy level
                        if (runif(1) > 0.5) {
                            ### Select the current allele as the other non-zero allele randomly
                            g[b] = 1 - sum(g[a:b])
                        }
                    }
                }
            }
            return(g)
        })
    }
	parallel::stopCluster(clusters)
    G = t(G)
    gc()
    ### Duplicate the vectors of loci coordinates by the number of $n_alleles-1$ alleles
    vec_chr = rep(vec_chr, each=n_alleles-1)
    vec_pos = rep(vec_pos, each=n_alleles-1)
    vec_allele = rep(paste0("allele_", 1:(n_alleles-1)), times=l)
    ### Output n entries x l loci matrix
    n_digits = length(unlist(strsplit(as.character(n), "")))
    rownames(G) = paste0("entry_", sprintf(paste0("%0", n_digits, "d"), 1:n))
    colnames(G) = paste(vec_chr, vec_pos, vec_allele, sep="-")
    if (verbose==TRUE) {
        print("Genotype simulation finished:")
        print(paste0("     -number of rows (number of entries) = ", n))
        print(paste0("     -number of columns (number of loci x number of alleles - 1) = ", l*(n_alleles-1)))
        print("Allele frequency distribution:")
        txtplot::txtdensity(G[sample(x=n, size=min(c(n, 1000))), sample(x=m, size=min(c(n, 1000)))])
        print("##############################################################################")
    }
    if (show_correlation_heatmap==TRUE) {
        idx = seq(from=2, to=ncol(G), by=(n_alleles-1))
        heatmap(cor(G[, idx]), Rowv=NA, Colv=NA)
    }
    return(G)
}

#' Simulate phenotype data from purely additive or additive + epistatic effects.
#' For purely additive effects, all the simulated QTL are used to compute the genetic value of each entry.
#' For additive + epistatic effects and if the number of effects per network required is less than the number of possible base-3 combinations of effects,
#' then the first 90% of the required effects per network corresponding to the first consecutive additive+epistatic effects are included and the remaining are sampled regularly for the remaining effect combinations.
#' This may result in less than the required effects per network, if the number of combinations are very large.
#' 
#' @param G: genotype matrix with $n$ rows, $l x (n_alleles-1)$ columns, and named rows and columns [Required]
#' @param n_alleles: number of alleles per locus [Required]
#' @param dist_effects: distribution of genetic effects which may be "norm" for Gaussian or "chi2" for Chi-squared
#' @param n_effects: number of additive genetic effects
#' @param purely_additive: Simulate only additive effects?
#' @param n_networks: number of networks with non-additive effects
#' @param n_effects_per_network: maximum number of non-additive effects per network
#' @param h2: heritability, i.e., variance due to genetic effects (additive or additive and non-additive) divided by total phenotype variance
#' @param pheno_reps: number of times the phenotypes will be recalculated where the error effects are resampled for each sample 
#' @param verbose: show simulation messages?
#' @returns
#' Y: phenotype matrix with $n$ rows, $pheno_reps$ columns, and named rows and columns
#' b: additive genetic effects, if purely_additive == TRUE, NULL otherwise
#' E: epistasis matrix, if purely_additive == FALSE, NULL otherwise
#' b_epi: epistasis effects, if purely_additive == FALSE, NULL otherwise
#' @examples
#' G = fn_simulate_genotypes()
#' list_Y_b_E_b_epi_1 = fn_simulate_phenotypes(G)
#' list_Y_b_E_b_epi_2 = fn_simulate_phenotypes(G=G, n_alleles=4, dist_effects="chi2", n_effects=25, h2=0.75, pheno_reps=5, verbose=TRUE)
#' list_Y_b_E_b_epi_3 = fn_simulate_phenotypes(G=G, n_effects=25, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, verbose=TRUE)
#' @export
fn_simulate_phenotypes = function(G, n_alleles=2, dist_effects=c("norm", "chi2")[1], n_effects=5, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.5, pheno_reps=1, verbose=FALSE) {
    # n_alleles=2
    # G = fn_simulate_genotypes(n_alleles=n_alleles)
    # dist_effects=c("norm", "chi2")[1]
    # n_effects=100
    # purely_additive=TRUE
    # n_networks=10
    # n_effects_per_network=50
    # h2=0.5
    # pheno_reps=1
    # verbose=TRUE
    ### Define the dimensions of the genotype data
    n = nrow(G)
    p = ncol(G)
    l = round(p / (n_alleles-1))
    ### Remove fixed loci
    idx_retain = c()
    for (j in 1:p) {
        if (var(G[,j]) >= 1e-7) {
            idx_retain = c(idx_retain, j)
        }
    }
    if (length(idx_retain) < p) {
        G = G[, idx_retain]
    }
    ### Define the probability distribution function from which the allele effects will be sampled from
    if (dist_effects=="norm") {
        fn = function(){rnorm(n=1, mean=0, sd=1)}
    } else if (dist_effects=="chi2") {
        fn = function(){rchisq(n=1, df=1)}
    } else {
        print("Please select: 'norm' or 'chi2'.")
    }
    ### Randomly select the QTL
    if (n_effects > l) {
        n_effects = l
    }
    idx_effect = sort(sample(x=1:l, size=n_effects))
    ### Define the allele effects across the QTL
    ### n_alleles - 1 alleles will be given effects, i.e. the minor allele per QTL wil have no effect for simplicity
    b = rep(0, times=p)
    for (j in idx_effect) {
        for (k in 0:(n_alleles-2)) {
            idx = (j-1)*(n_alleles-1) + 1 + k
            b[idx] = fn()
        }
    }
    ### Define the genetic effects
    if (purely_additive) {
        ### Define the breeding values
        g = G %*% b
        E = NULL
        b_epi = NULL
    } else {
        ### Define additive and epistatic or interaction effects
        ### The interaction or epistatic effects are simply the products of the genotypes x additive effects across multiple QTL,
        ### this means we can have 2-way, 3,-way, 4-way, etc... interaction effects.
        ### We control sparsity of these interaction effects (which grows exponential as the number of QTLs increases) by:
        ###     - dividing the QTLs into gene/allele networks, and by
        ###     - setting a maximum number of effects per network, i.e. n_effects_per_network.
        ### Randomise the sorting, and hence grouping of the QTLs
        idx_effect_randomised = sample(idx_effect, size=length(idx_effect), replace=FALSE)
        n_idx_per_network = floor(length(idx_effect_randomised) / n_networks)
        ### Iterate across gene networks
        if (verbose==TRUE) {
            print("Simulating additive and epistatic effects with:")
            print(paste0("     - ", n_effects, " QTLs distributed across"))
            print(paste0("     - ", n_networks, " gene networks with a maximum of"))
            print(paste0("     - ", n_effects_per_network, " effects per gene network."))
            pb = txtProgressBar(min=0, max=n_networks*n_effects_per_network, style=3)
        }
        for (i in 1:n_networks) {
            # i = 1
            ### Correct for the number of QTL included in each network when they don't divide evenly
            if ((i == n_networks) & ((n_idx_per_network * n_networks) < length(idx_effect_randomised))) {
                idx_randomised = idx_effect_randomised[((i-1)*n_idx_per_network+1):length(idx_effect_randomised)]
            } else {
                idx_randomised = idx_effect_randomised[((i-1)*n_idx_per_network+1):(i*n_idx_per_network)]
            }
            ### The total number of epistatic effects is base 3 because we have 3 columns initially: (1) additive effects of the previous allele, (2) additive effects of the current allele, and (3) the interaction between the two.
            m_max = 3^(length(idx_randomised)-1)
            ### If we are asking for more effects per network than is possible, then we set the number of effects per network to this maximum value.
            if (n_effects_per_network > m_max) {
                m_i = m_max
            } else {
                m_i = n_effects_per_network
            }
            ### We need to make sure the first 90% of the additive and epistatic effects we are asking for are completely covered so we don't run out of effects to select from subsequently.
            ### We then take regularly spaced samples of the additive and epistatic effects across the rest of the uncovered m_max effects amounting to 90% of the effects we're asking for.
            idx_of_epistatic_effects_included = floor(c(1:floor(0.9*m_i), seq(from=floor(0.9*m_i)+1, to=m_max, length=(m_i-floor(0.9*m_i)))))
            for (j in idx_randomised) {
                # j = idx_randomised[1]
                ### Calculate the additive effects per entry: allele genotype x allele effects.
                # x = G[, j, drop=FALSE] * b[j]
                x = G[, j, drop=FALSE]
                if (j==idx_randomised[1]) {
                    ### Initialise the matrix of additive and epistatic effects.
                    F = scale(x, scale=TRUE, center=TRUE)
                    # F = x
                    ### We use this counter to keep track of the expected number of additive and epistatic effects so that we can map it to the indexes of the effects we want to include.
                    counter = 1
                } else {
                    m = ncol(F)
                    for (k in 1:m) {
                        # k = 1
                        ### Extract the kth additive/epistatic effect we will be using multiplying with current additive effect, x
                        f = F[, k, drop=FALSE]
                        counter = counter + 1
                        if (sum(counter == idx_of_epistatic_effects_included) > 0) {
                            ### Append the current additive effect
                            F = cbind(F, x)
                        }
                        counter = counter + 1
                        if (sum(counter == idx_of_epistatic_effects_included) > 0) {
                            ### Append the interaction effect between the kth additive/epistatic effect and the current additive effect
                            rand_within = sample(x=c(1, 2, 3), size=2) ### none, 1/x, or ln(abs(x))
                            rand_across = sample(x=c(1, 2, 3), size=1) ### multiply, divide, power
                            if (rand_within[1]==1) {
                                a = f
                            } else if (rand_within[1]==2) {
                                a = sign(f)/(abs(f)+1e-7)
                            } else {
                                a = log(abs(f)+1e-7)
                            }
                            if (rand_within[2]==1) {
                                b = x
                            } else if (rand_within[2]==2) {
                                b = sign(x)/(abs(x)+1e-7)
                            } else {
                                b = log(abs(x)+1e-7)
                            }
                            if ((var(a) < 1e-7) | (var(b) < 1e-7)) {
                                next
                            }
                            a_min = min(a[is.infinite(a)==FALSE]); a_max = max(a[is.infinite(a)==FALSE])
                            a[(sign(a)==-1) & is.infinite(a)] = a_min; a[(sign(a)==+1) & is.infinite(a)] = a_max
                            b_min = min(b[is.infinite(b)==FALSE]); b_max = max(b[is.infinite(b)==FALSE])
                            b[(sign(b)==-1) & is.infinite(b)] = b_min; b[(sign(b)==+1) & is.infinite(b)] = b_max
                            a = scale(a, scale=TRUE, center=TRUE)
                            b = scale(b, scale=TRUE, center=TRUE)
                            if (rand_across==1) {
                                e = a*b
                            } else if (rand_across==2) {
                                e = sign(b)*a/(abs(b)+1e-7)
                            } else {
                                e = abs(a)^b
                            }
                            e_min = min(e[is.infinite(e)==FALSE]); e_max = max(e[is.infinite(e)==FALSE])
                            e[(sign(e)==-1) & is.infinite(e)] = e_min; e[(sign(e)==+1) & is.infinite(e)] = e_max
                            # e = scale(e, scale=TRUE, center=TRUE)
                            F = cbind(F, e)
                            colnames(F)[ncol(F)] = paste0(colnames(f), "-*-", colnames(x))
                        }
                    }
                }
                if (verbose==TRUE) {
                    ### Update progress bar
                    if (i==1) {
                        setTxtProgressBar(pb, ncol(F))
                    } else {
                        setTxtProgressBar(pb, ncol(E)+ncol(F))
                    }
                }
                ### Update counter to the expected number of columns if there were no selection for only a subset of m_max
                counter =  3^(which(idx_randomised==j)-1)
            }
            if (verbose==TRUE) {
                ### Update the progress bar with the target number of columns (we may not have reached this target when m_max is very large)
                setTxtProgressBar(pb, i*m_i)
            }
            ### Build the additive and epistasis effects matrix
            if (i==1) {
                E = F
            } else {
                E = cbind(E, F)
            }
        }
        if (verbose==TRUE) {
            setTxtProgressBar(pb, n_networks*n_effects_per_network)
            close(pb)
        }
        ### Define the epistatic effects
        b_epi = c()
        for (j in 1:ncol(E)) {
            b_epi = c(b_epi, fn())
        }
        ### Define the genotype effects as the sum across additive and epistatic effects without normalisation to avoid wiping out the non-linear effects
        g = E %*% b_epi
        b = NULL
    }
    ### Define the residual variance but first we standard normalise the genotype values for simplicity
    g = scale(g, scale=TRUE, center=TRUE)
    vg = var(g)
    ve = 1/h2 - 1
    ### Define the phenotypes for each replicate, where each replicate is associated with independently sampled residual effects
    for (i in 1:pheno_reps) {
        residual = rnorm(n, mean=0, sd=sqrt(ve))
        y_i = g + residual
        if (i == 1) {
            Y = matrix(y_i, ncol=1)
        } else {
            Y = cbind(Y, y_i)
        }
    }
    if (verbose==TRUE) {
        txtplot::txtdensity(Y, xlab="Phenotype Values", ylab="Density")
    }
    ### Output the phenotype values and genotype values
    rownames(Y) = rownames(G)
    colnames(Y) = paste0("rep_", 1:pheno_reps)
    return(list(Y=Y, b=b, E=E, b_epi=b_epi))
}

#' Simulate genotype-by-environment interactions
#' Arguments are the same as fn_simulate_phenotypes excluding pheno_reps which is replaced with n_reps which refers to the number of replicates of observations per genotype per environment in addition to env_factor_levels, and env_factor_effects_sd.
#' 
#' @param G: genotype matrix with $n$ rows, $l x (n_alleles-1)$ columns, and named rows and columns [Required]
#' @param n_alleles: number of alleles per locus [Required]
#' @param dist_effects: distribution of genetic effects which may be "norm" for Gaussian or "chi2" for Chi-squared
#' @param n_effects: number of additive genetic effects
#' @param purely_additive: Simulate only additive effects?
#' @param n_networks: number of networks with non-additive effects
#' @param n_effects_per_network: maximum number of non-additive effects per network
#' @param h2: heritability, i.e., variance due to genetic effects (additive or additive and non-additive) divided by total phenotype variance
#' @param env_factor_levels: vector of environmental factor levels - think of this as the number of distinct classes per environmental factor like precipitation, temperature, and solar radiation regimes
#' @param env_factor_effects_sd: vector or a single value defining the standard deviation of the normal distribution centred at 0 from which the environmental effects affecting the genotype effect will be sampled from
#' @param n_reps: number of replicates of observations per genotype per environment
#' @param verbose: show simulation messages?
#' @returns 
#' df: a data.frame with the phenotype in a single column, replication ID, genotype ID, environment ID, followed by the environmental factor levels with one column per factor
#' CORR: correlation matrix between the environments
#' @examples
#' G = fn_simulate_genotypes()
#' list_df_CORR_1 = fn_simulate_gxe(G)
#' list_df_CORR_2 = fn_simulate_gxe(G=G, env_factor_levels=c(5, 3), env_factor_effects_sd=0.2, n_reps=5, verbose=TRUE)
#' list_df_CORR_3 = fn_simulate_gxe(G=G, n_effects=50, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, h2=0.5, env_factor_levels=c(5, 3), env_factor_effects_sd=0.2, n_reps=5, verbose=TRUE)
#' @export
fn_simulate_gxe = function(G, n_alleles=2, dist_effects=c("norm", "chi2")[1], n_effects=5, purely_additive=TRUE, n_networks=1, n_effects_per_network=50, h2=0.5, env_factor_levels=c(2, 3, 2), env_factor_effects_sd=c(0.1, 1.0, 0.01), frac_additional_QTL_per_env=0.01, n_reps=3, verbose=FALSE) {
    # G = fn_simulate_genotypes()
    # n_alleles = 2
    # dist_effects=c("norm", "chi2")[1]
    # n_effects=5
    # purely_additive=FALSE
    # # purely_additive=TRUE
    # n_networks=10
    # n_effects_per_network=50
    # h2=0.5
    # env_factor_levels=c(2, 3, 2)
    # env_factor_effects_sd=c(0.1, 1.0, 0.01)
    # frac_additional_QTL_per_env=0.01
    # n_reps=3
    # verbose=FALSE
    n = nrow(G)
    l = ncol(G)
    ### Simulate genetic effects
    list_Y_b_E_b_epi = fn_simulate_phenotypes(G=G, n_alleles=n_alleles, dist_effects=dist_effects, n_effects=n_effects, purely_additive=purely_additive, n_networks=n_networks, n_effects_per_network=n_effects_per_network, h2=h2, pheno_reps=1, verbose=verbose)
    if (purely_additive==FALSE) {
        G = list_Y_b_E_b_epi$E
    }
    ### Instantiate output data.frame components
    y = c()
    rep = c()
    gen = c()
    env = c()
    n_env_factors = length(env_factor_levels)
    list_env_factors = eval(parse(text=paste0("list(", paste(paste0("env_factor_", 1:n_env_factors, "=c()"), collapse=", "), ")")))
    ### Define the standard deviation of the normal distribution centred at 0 from which the environmental effects affecting the genotype effect will be sampled from
    if (length(env_factor_effects_sd) < n_env_factors) {
        env_factor_effects_sd = rep(env_factor_effects_sd, times=n_env_factors)
    }
    ### Define new phenotypes per environment
    df_factor_level_combinations = eval(parse(text=paste0("expand.grid(", paste(paste0("env_factor_", 1:n_env_factors, "=c(1:", env_factor_levels, ")"), collapse=", "), ")")))
    for (i in 1:nrow(df_factor_level_combinations)) {
        # i = 1
        for (r in 1:n_reps) {
            if (purely_additive==TRUE) {
                effects = list_Y_b_E_b_epi$b
            } else {
                effects = list_Y_b_E_b_epi$b_epi
            }
            for (j in 1:ncol(df_factor_level_combinations)) {
                # j = 1
                effects = effects * df_factor_level_combinations[i, j] * rnorm(n=length(effects), mean=0, sd=env_factor_effects_sd[j])
            }
            ### Define the residual variance
            n_effects = sum(effects != 0.0)
            n_additional_effects = ceiling(frac_additional_QTL_per_env*n_effects)
            if ((n_effects + n_additional_effects) < length(effects)) {
                ### Define the probability distribution function from which the allele effects will be sampled from
                if (dist_effects=="norm") {
                    fn = function(){rnorm(n=1, mean=0, sd=1)}
                } else if (dist_effects=="chi2") {
                    fn = function(){rchisq(n=1, df=1)}
                } else {
                    print("Please select: 'norm' or 'chi2'.")
                }
                idx_additional_effects = sample(x=which(effects == 0.0), size=n_additional_effects)
                for (idx_new_effects in idx_additional_effects) {
                    effects[idx_new_effects] = fn()
                }
            }
            g = G %*% effects
            vg = var(g)
            ve = vg * (1/h2 - 1)
            y_new = scale(g + rnorm(n=n, mean=0, sd=ve), scale=TRUE, center=TRUE) ### Standard normalise the phenotypes
            y = c(y, y_new)
            rep = c(rep, rep(r, times=n))
            gen = c(gen, rownames(G))
            env = c(env, rep(paste("ENV", paste(df_factor_level_combinations[i, ], collapse="-"), sep="-"), times=n))
            for (j in 1:ncol(df_factor_level_combinations)) {
                list_env_factors[[j]] = c(list_env_factors[[j]], rep(df_factor_level_combinations[i, j], times=n))
            }
        }
    }
    ### Build the output data.frame
    df = data.frame(y, rep=as.factor(rep), gen=as.factor(gen), env=as.factor(env), as.data.frame(list_env_factors))
    ### Estimate correlation between environments
    CORR = matrix(NA, nrow=nlevels(df$env), ncol=nlevels(df$env))
    for (i in 1:nlevels(df$env)) {
        env_1 = levels(df$env)[i]
        for (j in 1:nlevels(df$env)) {
            env_2 = levels(df$env)[j]
            y_1 = df$y[df$env==env_1]
            y_2 = df$y[df$env==env_2]
            CORR[i, j] = cor(y_1, y_2)
            if (verbose & (i != j) & ((CORR[i, j] >= 0.5) | (CORR[i, j] <= -0.5))) {
                txtplot::txtplot(y_1, y_2, xlab=env_1, ylab=env_2)
                print(CORR[i, j])
            }
        }
    }
    rownames(CORR) = levels(df$env)
    colnames(CORR) = levels(df$env)
    if (verbose) {
        txtplot::txtdensity(df$y)
    }
    ### Output
    return(list(df=df, CORR=CORR))
}

#' Simulate sparse mixed-model data
#' 
#' @param G: genotype matrix with $n$ rows, $l x (n_alleles-1)$ columns, and named rows and columns [Required]
#' @param df: a data.frame with the phenotype in a single column, replication ID, genotype ID, environment ID, followed by the environmental factor levels with one column per factor
#' @param frac_gen_missing: fraction of the genotypes to be set as missing
#' @returns 
#' y_complete: vector of phenotype data without missing data
#' y: vector of phenotype data without simulated missing data
#' X: incidence matrix for the fixed effects including the intercept and environment/s
#' Z: incidence matrix for the random effects, i.e. the genotype IDs
#' D: variance-covariance matrix of the random effects defined as the Gram matrix of the genotype information, i.e. let $G$ be the numeric $n \times l$ genotype matrix which may represent the presence or absence of a genome marker or allele frequencies, then $D = GG^T \over l$.
#' @examples
#' G = fn_simulate_genotypes()
#' df_gxe = fn_simulate_gxe(G=G, n_effects=50, purely_additive=FALSE, n_networks=10, n_effects_per_network=50, h2=0.5, env_factor_levels=c(5, 3), env_factor_effects_sd=0.2, n_reps=5, verbose=TRUE)$df
#' list_y_complete_y_X_Z_D = fn_simulate_sparse_mixed_model_data(G=G, df_gxe=df_gxe, frac_gen_missing=0.25)
#' @export
fn_simulate_sparse_mixed_model_data = function(G, df_gxe, frac_gen_missing=0.1) {
    # G = fn_simulate_genotypes()
    # df_gxe = fn_simulate_gxe(G=G, env_factor_levels=c(1))$df
    # frac_gen_missing = 0.1
    ### Phenotype vector
    n = nrow(df_gxe)
    y = matrix(df_gxe$y, nrow=n, ncol=1)
    rownames(y) = paste0(df_gxe$gen, "-env_", gsub("-", "", df_gxe$env), "-rep_", df_gxe$rep)
    ### Fixed effects matrix (intercept and environment IDs)
    vec_env = sort(unique(as.character(df_gxe$env)))
    p = length(vec_env)
    X = matrix(1, nrow=n, ncol=1)
    if (p > 1) {
        for (j in 1:p) {
            # j = 1
            x = as.numeric(df_gxe$env == vec_env[j])
            X = cbind(X, x)
        }
        rownames(X) = rownames(y)
        colnames(X) = c("intercept", vec_env)
    } else {
        rownames(X) = rownames(y)
        colnames(X) = "intercept"
    }
    ### Random effects matrix (genotype IDs)
    vec_gen = sort(unique(as.character(df_gxe$gen)))
    q = length(vec_gen)
    Z = matrix(0, nrow=n, ncol=q)
    for (i in 1:n) {
        # i = 1
        j = which(df_gxe$gen[i] == vec_gen)
        Z[i, j] = 1
    }
    rownames(Z) = rownames(y)
    colnames(Z) = vec_gen
    ### Variance-covariance matrix of the random effects   
    D = G%*%t(G) / ncol(G)
    rownames(D) = colnames(Z)
    colnames(D) = colnames(Z)
    ### Simulate missing genotypes
    n_gen_missing = floor(frac_gen_missing * q)
    vec_gen_missing = sort(sample(vec_gen, size=n_gen_missing, replace=FALSE))
    vec_idx_gen_missing = df_gxe$gen %in% vec_gen_missing
    y_sparse = y
    y_sparse[vec_idx_gen_missing] = NA
    return(
        list(y_complete=y, y=y_sparse, X=X, Z=Z, D=D)
    )
}
