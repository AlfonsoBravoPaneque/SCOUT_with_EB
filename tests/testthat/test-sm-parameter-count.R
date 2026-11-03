test_that('SM does not count fixed tau for BM1, OU1, EB or OUM', {
    fit <- function(model, skip_tau, n_theta = 1L) {
        paras <- list(tau = 0, sigma = 1, theta = rep(2, n_theta))
        if (model == 'BM1') paras$alpha <- 1e-10
        if (model %in% c('OU1', 'OUM')) paras$alpha <- 0.5
        if (model == 'EB') paras$beta <- -0.5
        list(settings = c('gene', model), paras = paras,
             final = data.frame(iter = 1L, ll_total = -10, tau = 0),
             defaults = list(model = model, skipTau = skip_tau))
    }

    for (model in c('BM1', 'OU1', 'EB', 'OUM')) {
        n_theta <- if (model == 'OUM') 2L else 1L
        em <- fit(model, FALSE, n_theta)
        sm <- fit(model, TRUE, n_theta)
        em_k <- if (model == 'BM1') 3L else if (model == 'OUM') 5L else 4L

        em_history <- SCOUT:::extract_history_grid_search(list(em))
        sm_history <- SCOUT:::extract_history_grid_search(list(sm))
        expect_equal(em_history$param.count, em_k, info = model)
        expect_equal(sm_history$param.count, em_k - 1L, info = model)

        em_aic <- SCOUT:::fit_aic(em, 25L)
        sm_aic <- SCOUT:::fit_aic(sm, 25L)
        expect_equal(em_aic$k, em_k, info = model)
        expect_equal(sm_aic$k, em_k - 1L, info = model)
        expect_equal(sm_aic$aic, -2 * sm$final$ll_total + 2 * (em_k - 1L))
        expect_equal(sm_aic$aicc, -2 * sm$final$ll_total +
                         2 * (em_k - 1L) * 25 / (25 - (em_k - 1L) - 1))
    }
})
