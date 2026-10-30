test_that('EB covariance matches transformed tree and tends to BM at beta = 0', {
    for (nwk in c('((a:1,b:1):1,(c:1,d:1):1);',
                  '((a:0.5,b:1):1,(c:1,d:2):0.5);')) {
        phy <- ape::read.tree(text = nwk)
        depths <- ape::node.depth.edgelength(phy)
        info <- list(n_leaves = 4L, leaf_dists = depths[1:4],
                     shared_lengths = matrix(depths[ape::mrca(phy)], 4, 4))
        for (beta_H in c(-5, -1, -0.01, -1e-9, 0)) {
            beta <- beta_H / max(info$leaf_dists)
            actual <- SCOUT:::compute_VCV_EB(info, beta, 1.3)
            tr <- SCOUT:::transformPhy_EB(phy, beta, 1.3)
            expect_equal(actual, unname(ape::vcv(tr$tree)), tolerance = 1e-8)
            expect_equal(tr$tree$root.edge, 0)
            expect_equal(unname(tr$diagWeight), rep(1, 4))
        }
        expect_equal(SCOUT:::compute_VCV_EB(info, 0, 1.3),
                     unname(1.3 * ape::vcv(phy)))
        expect_error(SCOUT:::compute_VCV_EB(info, 0.1, 1.3), 'beta <= 0', fixed = TRUE)
        expect_error(SCOUT:::transformPhy_EB(phy, 0.1, 1.3), 'beta <= 0', fixed = TRUE)
        scaled <- info
        scaled$leaf_dists <- 100 * info$leaf_dists
        scaled$shared_lengths <- 100 * info$shared_lengths
        expect_equal(SCOUT:::compute_VCV_EB(info, -0.3, 1.3),
                     SCOUT:::compute_VCV_EB(scaled, -0.3 / 100, 1.3 / 100))
    }
})

test_that('EB is a single-regime model with fixed root regardless of user default', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    meta <- data.frame(species = phy$tip.label, gene = c(2, 3, 4, 5))
    idata <- SCOUT::formatSCOUT(phy, meta, tempdir(), regimes = c('BM1', 'OU1', 'EB'),
                                anc_infer = 'skip', normalize = FALSE)
    expect_equal(idata$gene_cols, 'gene')
    for (fixed in c(FALSE, TRUE)) {
        ti <- SCOUT:::preprocessTree(idata, 'EB', root.fixed = fixed)
        expect_true(ti$defaults$add.root)
        expect_equal(ti$defaults$model, 'EB')
        gi <- SCOUT:::preprocessGene(idata, 'gene', ti$tree_const, ti$defaults$root.state,
                                    ti$defaults$add.root, 'EB', FALSE)
        expect_lte(gi$param_init$beta, 0)
        expect_null(gi$param_init$alpha)
        expect_length(gi$param_init$theta, 1)
    }
})

test_that('EB near-BM flag uses beta times tree height without forcing a BM call', {
    rows <- data.frame(gene_name = 'g', regime = c('BM1', 'OU1', 'EB'),
                       model = c('BM1', 'OU1', 'EB'), ll_total = c(-100, -100, -95),
                       param.count = c(3, 4, 4), ntips = 100,
                       beta = c(NA, NA, -0.0001), beta_H = c(NA, NA, -0.01),
                       converge = 'converged', iter = 8L)
    res <- SCOUT::scout_model_selection(rows, hybrid = TRUE)
    expect_equal(res$selected_regime, 'EB')
    expect_true(res$beta_near_bm)
    expect_equal(res$lrt_status, 'not eligible')
    expect_equal(res$beta_H, -0.01)
    rows$beta_H[3] <- -0.011
    expect_false(SCOUT::scout_model_selection(rows, hybrid = TRUE)$beta_near_bm)
})

test_that('hybrid LRT excludes EB from OU multi-regime tests and fallback', {
    rows <- data.frame(regime = c('BM1', 'OU1', 'EB', 'OU2'),
                       ll_total = c(-120, -110, -109.9, -108.8),
                       param.count = c(3, 4, 4, 5), ntips = 100)
    res <- SCOUT::select_model_hybrid(rows, ic = 'AIC')
    expect_equal(res$ic_winner, 'OU2')
    expect_equal(res$n_multi, 1L)
    expect_equal(res$best_simple, 'EB')
    expect_equal(res$call, 'EB')
    expect_equal(res$lrt_alt, 'OU2')
    expect_equal(res$lrt_null, 'OU1')
})

test_that('EB M-step and E-step use the same covariance and count beta', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    meta <- data.frame(species = phy$tip.label, gene = c(2, 3, 4, 5))
    idata <- SCOUT::formatSCOUT(phy, meta, tempdir(), regimes = 'EB',
                                anc_infer = 'skip', normalize = FALSE)
    ti <- SCOUT:::preprocessTree(idata, 'EB', root.fixed = FALSE)
    gi <- SCOUT:::preprocessGene(idata, 'gene', ti$tree_const, ti$defaults$root.state,
                                ti$defaults$add.root, 'EB', FALSE)
    d <- ti$defaults
    d$verbose <- FALSE
    d$tau_prior_mean <- 0.2
    d$tau_prior_sd <- 0.1
    er <- SCOUT:::e_step(gi$X, ti$tree_const$tree, ti$tree_const$edges,
                         gi$param_init, d, diagnose = FALSE)
    expect_length(er$E_Z, 4)
    pred <- SCOUT:::predict_masked_latent(gi$param_init, d, seq_len(4),
                                           gi$X, kept_idx = seq_len(4))
    expect_equal(unname(as.vector(er$E_Z)), unname(pred$mean), tolerance = 1e-8)
    fit <- SCOUT:::m_step_tipfog(gi$X, gi$X, gi$param_init,
                                 ti$tree_const$tree, ti$tree_const$edges, d, 0.2, 0.2)
    expect_lte(fit$beta, 0)
    V <- SCOUT:::compute_VCV_EB(d$parsed_alt_tree, fit$beta, fit$sigma)
    diag(V) <- diag(V) + fit$tau^2
    expect_equal(fit$ll_total,
                 as.numeric(-0.5 * (4 * log(2 * pi) + as.numeric(determinant(V, logarithm = TRUE)$modulus) +
                                    t(gi$X - fit$theta) %*% solve(V, gi$X - fit$theta))),
                 tolerance = 1e-7)
    out <- SCOUT:::run_em(ti$tree_const$tree, gi$X, gi$param_init,
                          ti$tree_const$edges, d, max_iter = 1, skipE = TRUE,
                          save = FALSE, diagnose = FALSE, verbose = FALSE)
    expect_equal(length(unlist(out$paras)), 4)
    expect_equal(out$final$beta_H,
                 out$paras$beta * max(d$parsed_alt_tree$leaf_dists))
})

test_that('EB passes through runSCOUT parameter extraction and IC accounting', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    meta <- data.frame(species = phy$tip.label, gene = c(2, 3, 4, 5),
                       gene2 = c(3, 2, 5, 4))
    idata <- SCOUT::formatSCOUT(phy, meta, tempdir(), regimes = 'EB',
                                anc_infer = 'skip', normalize = FALSE)
    fits <- SCOUT::runSCOUT(idata, lambda1 = 0.2, lambda2 = 0.2,
                            fixed.root = FALSE, M_only = TRUE, cores = 1, verbose = FALSE)
    hist <- SCOUT:::extract_history_grid_search(fits)
    expect_true(all(hist$param.count == 4))
    expect_true(all(hist$model == 'EB'))
    pars <- SCOUT::extract_parameters(fits)
    expect_equal(length(pars$EB$gene$theta_1), 1)
    expect_true(is.finite(pars$EB$gene$beta))
    expect_equal(hist$beta_H,
                 hist$beta * max(fits[[1]]$defaults$parsed_alt_tree$leaf_dists))
    idata$gene_cols <- 'gene'
    single <- SCOUT::runSCOUT(idata, lambda1 = 0.2, lambda2 = 0.2,
                              fixed.root = FALSE, M_only = TRUE, cores = 1, verbose = FALSE)
    expect_length(single, 1)
})
