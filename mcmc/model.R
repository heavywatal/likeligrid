#########1#########2#########3#########4#########5#########6#########7#########
# Exhaustive loops
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

.results = .grid %>>%
#    head() %>>%
    purrr::by_row(~calc_loglik(.genotype_cancer, flatten_dbl(.x)), .to='loglik') %>>%
    tidyr::unnest() %>>%
    dplyr::arrange(desc(loglik)) %>>%
    (?.)

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

#########1#########2#########3#########4#########5#########6#########7#########
# Visualize

rescale = function(x, from=0, to=1) {
    .range = range(x)
    (x - .range[1]) / (.range[2] - .range[1])
}

.forplot = .results %>>%
    dplyr::filter(loglik > 1.1 * max(loglik)) %>>%
    dplyr::mutate(
        x= rescale(loglik),
        color= rgb(x, 1 - x, 1 - x)) %>>% (?.)

#########1#########2#########3#########4

library(rgl)

if (rgl.cur()) {rgl.close()}
rgl::open3d(windowRect=c(0, 0, 600, 600))
.forplot %>>%
    {rgl::spheres3d(x=.$TGF, y=.$Damage, z=.$Cycle,
       color=.$color, radius=0.03)}
rgl::axes3d(box=TRUE, xat=seq(0, 10, 0.05))
rgl::box3d(xat=seq(0, 100, 0.1), yat=c(0, 2), zat=c(0, 1))
#rgl::title3d('', '', 'x', 'y', 'z')
rgl::writeWebGL('.', 'loglik_landscape.html', snapshot=FALSE, width=600, height=600)

#########1#########2#########3#########4

library(plotly)

.forplot %>>%
    plot_ly(x=~TGF, y=~Damage, z=~Cycle, type='scatter3d', mode='markers', marker=list(color=~color, size=24))

#########1#########2#########3#########4#########5#########6#########7#########
# deprecated

unnest_pathway = function(.data) {
    dplyr::mutate(.data, pathway= str_split(pathway, ':')) %>>%
    tidyr::unnest() %>>%
    dplyr::arrange(sample, pathway, desc(mut_score))
}

weight_by_gene1 = function(.data, .weight=init_genes) {
    dplyr::mutate(.data, mut_score= mut_score * .weight[gene])
}

apply_gene_model = function(.data, .func) {
    dplyr::group_by(.data, sample, pathway) %>>%
    dplyr::summarise(mut_score= .func(mut_score))
}

weight_by_pathway1 = function(.data, .weight=init_pathways) {
    dplyr::mutate(.data, mut_score= mut_score * .weight[pathway])
}

calc_loglik = function(.data, sd_mut=10, is_cancer=TRUE) {
    dplyr::summarise(.data, mut_score= sum(mut_score), n_mut=n()) %>>%
    dplyr::mutate(loglik= pnorm(100, mut_score, sd_mut * n_mut, lower.tail=!is_cancer, log.p=TRUE))
}

maf_pathways = maf_weighted %>>%
    unnest_pathway() %>>%
    weight_by_gene1() %>>%
    apply_gene_model(max) %>>% (?.)

model1 = function(.data, .weight, .sd) {
    weight_by_pathway1(.data, .weight) %>>%
    calc_loglik(.sd) %>>%
    (sum(.$loglik))
}
model1(maf_pathways, init_pathways * 10, 10)


#########1#########2#########3#########4#########5#########6#########7#########
# Huge exhaustive loops

.by = 0.1
.v = seq(.by, 0.5, .by)
(.param_names = names(tbl_cancer))
(num_params = length(.param_names))

.by = 0.05
.v = seq(.by, 0.5, .by)
.grid = tidyr::crossing(V1=.v, V2=.v) %>>%
tidyr::expand(nesting(V1, V2), V3=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3), V4=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4), V5=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5), V6=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5, V6), V7=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5, V6, V7), V8=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5, V6, V7, V8), V9=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5, V6, V7, V8, V9), V10=V1) %>>%
dplyr::filter(rowSums(.) <= 1) %>>%
tidyr::expand(nesting(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10), V11=V1) %>>%
dplyr::filter(rowSums(.) == 1) %>>% (?.)

names(.grid) = .param_names

.sigma = 0.03
.sigmas = rep(.sigma, length(.param_names))
.threshold = 0.3

calc_loglik_cancer = function(.data, params, sigmas, threshold) {
    .genotypes = as.matrix(.data[,seq_along(params)])
    detail_calc_loglik(.genotypes, TRUE, params, sigmas, threshold) %>>%
    (bind_cols(.data, .)) %>>%
    (sum(.$loglik))
}
calc_loglik_cancer(tbl_cancer, flatten_dbl(.grid[3,]), .sigmas, .threshold)

.results_cancer = .grid %>>%
    purrr::by_row(.to='loglik', ~{
      calc_loglik_cancer(tbl_cancer, flatten_dbl(.x), .sigmas, .threshold)
    }) %>>%
    dplyr::mutate(loglik= flatten_dbl(loglik)) %>>%
    (?.)

.results03 = .results_cancer %>>% arrange(desc(loglik))

seq(0.3, 0.6, 0.1) %>>%
purrr::map(~{

})
