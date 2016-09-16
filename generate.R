#########1#########2#########3#########4#########5#########6#########7#########
## Make dummy data

n_obs = 4000L
true_threshold = 0.5
true_epsilon = 0.1
true_coefs = c(RAS=0.15, TGF=0.30, Notch=0.05, Cycle=0.25, Damage=0.15)
stopifnot(sum(true_coefs, true_epsilon) == 1.0)

genotypes = rbernoulli(n_obs * length(true_coefs), 0.5) %>>%
    as.integer() %>>%
    matrix(nrow=n_obs) %>>%
    (? head(.))

.data = genotypes %>>% {
        .mu = . %*% true_coefs + true_epsilon
        tibble::as_tibble(.) %>>%
        setNames(names(true_coefs)) %>>%
        dplyr::mutate(true_mu=.mu, true_sigma=sqrt(true_mu))
    } %>>%
    dplyr::mutate(
        true_score= rnorm(nrow(.), true_mu, true_sigma),
        cancer= true_score >= true_threshold) %>>%
    (?.) %>>% (? dplyr::count(., cancer))

.geno_matrices = .data %>>%
    dplyr::select(-matches('^true_')) %>>%
    tidyr::nest(-cancer) %>>%
    dplyr::mutate(data=purrr::map(data, as.matrix)) %>>% (?.)

.genotype_cancer = dplyr::filter(.geno_matrices, cancer)$data[[1]] %>>%
    (?as_tibble(.))

.genotype_cancer %>>% wtl::write_df('genotype_isovar_cancer.tsv')

#########1#########2#########3#########4#########5#########6#########7#########
# Use only cancer data

.epsilon = 0.1
.threshold = 0.5
.grid = 19
.v = seq(0.0, 1.0 - .epsilon, length.out=.grid)
num_params = ncol(.genotype_cancer)

.grid = rep(.v, num_params - 1) %>>%
    matrix(ncol=num_params - 1) %>>%
    tibble::as_tibble() %>>%
    expand.grid() %>>%
    tibble::as_tibble() %>>%
    dplyr::mutate(rest= 1.0 - .epsilon - rowSums(.)) %>>%
    dplyr::filter(rest >= 0) %>>%
    setNames(names(true_coefs)) %>>%
    (?.)

(.params = .grid %>>% sample_n(1) %>>% unlist())

calc_loglik = function(genotype_mat, params, cancer=TRUE) {
    dplyr::tibble(
      mu= genotype_mat %*% params %>>% as.vector() + .epsilon,
      sigma= sqrt(mu),
      loglik= pnorm(.threshold, mu, sigma, lower.tail=!cancer, log.p=TRUE)
    ) %>>% (sum(.$loglik))
}

calc_loglik(.genotype_cancer, .params)

.grid %>>%
#    head() %>>%
    purrr::by_row(~calc_loglik(.genotype_cancer, flatten_dbl(.x)), .to='loglik') %>>%
    tidyr::unnest() %>>%
    dplyr::arrange(desc(loglik))

#########1#########2#########3#########4
# Use non-cancer too

.grid %>>%
    purrr::by_row(~{
        .params = flatten_dbl(.x)
        purrr::by_row(.geno_matrices, ~{
            calc_loglik(.$data[[1]], .params, .$cancer)
        }, .labels=FALSE, .collate='rows')$.out %>>% sum()
    }, .to='loglik') %>>%
    tidyr::unnest() %>>%
    dplyr::arrange(desc(loglik))
