test_that('EB latent draws follow the fitting covariance, including the BM limit', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    for (beta in c(-1, 0)) {
        set.seed(42)
        draws <- SCOUT:::simulate_EB_EVFs(phy, beta, 0.5, 5, 4000)
        expected <- unname(ape::vcv(SCOUT:::transformPhy_EB(phy, beta, 0.5)$tree))
        expect_equal(unname(stats::cov(t(draws))), expected, tolerance = 0.06)
        expect_equal(rowMeans(draws), rep(5, 4), tolerance = 0.02)
        if (beta == 0) expect_equal(expected, unname(0.5 * ape::vcv(phy)))
    }
})

test_that('EB simulation retains OU output and produces aligned counts and latent values', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    phy$states <- c('A', 'A', 'B', 'B')
    phy$node.label <- c('A', 'A', 'B')
    od <- tempfile()
    dir.create(od)
    on.exit(unlink(od, recursive = TRUE))

    set.seed(1)
    sim <- simulate_test_data(ngenes = 4, ncells = 4, tree = phy, outdir = od,
                              out_prefix = 'eb', a = 1, s = 0.5, t0 = 5,
                              theta_step = 1, nevfs = 3, sim_lin = FALSE,
                              beta = c(-1, 0))
    expect_named(sim, c('scLT', 'simulated_OU', 'simulated_EB'))
    expect_length(sim$simulated_EB$counts, 2)
    expect_length(sim$simulated_EB$evfs, 2)
    for (cnt in sim$simulated_EB$counts) {
        expect_equal(cnt$species, phy$tip.label)
        expect_named(cnt, c(paste0('EB_', 1:4), 'species'))
        expect_true(all(as.matrix(cnt[, 1:4]) >= 0))
        expect_true(all(as.matrix(cnt[, 1:4]) == floor(as.matrix(cnt[, 1:4]))))
    }
    for (evf in sim$simulated_EB$evfs) {
        expect_equal(evf$species, phy$tip.label)
        expect_equal(ncol(evf), 10)
        expect_true(all(is.finite(as.matrix(evf[, 1:9]))))
    }
    expect_true(file.exists(file.path(od, 'eb_beta_-1_sigma_0.5_EB_counts.csv')))
    expect_true(file.exists(file.path(od, 'eb_beta_0_sigma_0.5_EB_evfs.csv')))

    set.seed(1)
    legacy <- simulate_test_data(ngenes = 4, ncells = 4, tree = phy, outdir = od,
                                 out_prefix = 'legacy', a = 1, s = 0.5, t0 = 5,
                                 theta_step = 1, nevfs = 3, sim_lin = FALSE)
    expect_named(legacy, c('scLT', 'simulated_OU'))
    expect_equal(legacy$simulated_OU, sim$simulated_OU)
})

test_that('EB simulation validates parameters and branch lengths', {
    phy <- ape::read.tree(text = '((a:1,b:1):1,(c:1,d:1):1);')
    for (beta in list(1, NA_real_, Inf, numeric())) {
        expect_error(simulate_test_data(2, 4, phy, tempdir(), 'invalid', 1, 1, 5, 1,
                                        beta = beta), 'finite values <= 0')
    }
    phy$edge.length <- NULL
    expect_error(SCOUT:::simulate_EB_genes(phy, -1, 1, 5, 2, 2, 4), 'branch lengths')
})
