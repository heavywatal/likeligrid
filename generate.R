#########1#########2#########3#########4#########5#########6#########7#########
## Functions

make_genotypes = function(n_obs, param_names, mu) {
    rbernoulli(n_obs * length(param_names), mu) %>>%
    as.integer() %>>%
    matrix(nrow=n_obs) %>>%
    {colnames(.) = param_names; .}
}

add_interaction_term = function(.data, .indices) {
    .names = names(.data)[.indices]
    .eq = paste(.names, collapse='*')
    names(.eq) = paste(.names, collapse=':')
    dplyr::mutate_(.data, .dots=.eq)
}
# as_tibble(genotypes) %>>% add_interaction_term(c(1, 2))
# as_tibble(genotypes) %>>% add_interaction_term(c(1, 2, 3))

make_genotypes2 = function(genotypes) {
    .dimensions = ncol(genotypes)
    .combn = combn(.dimensions, 2)
    .genotypes2 = as_tibble(genotypes)
    for (i in seq_len(ncol(.combn))) {
        .genotypes2 = add_interaction_term(.genotypes2, .combn[,i])
    }
    .genotypes2
}

make_phenotypes = function(genotypes, coefs, sd_coef, threshold, epsilon=0.1) {
    genotypes %>>% {
        .mu = . %*% coefs + epsilon
        tibble::as_tibble(.) %>>%
        dplyr::mutate(true_mu=.mu)
    } %>>% dplyr::mutate(
        true_sigma=sqrt(true_mu) * sd_coef,
        true_score= rnorm(nrow(.), true_mu, true_sigma),
        cancer= true_score >= threshold)
}

#########1#########2#########3#########4#########5#########6#########7#########
## Make dummy data

n_obs = 16000L
mu = 0.05
true_threshold = 0.25
true_epsilon = 0.1
sd_coef = 0.2
true_coefs = c(0.15, 0.30, 0.05, 0.25, 0.15)
names(true_coefs) = head(names(genotype_pereira), length(true_coefs))
stopifnot(sum(true_coefs, true_epsilon) == 1.0)

genotypes = make_genotypes(n_obs, names(true_coefs), mu) %>>%
    (? head(.))

.data = make_phenotypes(genotypes, true_coefs, sd_coef, true_threshold, true_epsilon) %>>%
    (?.) %>>% (? dplyr::count(., cancer))

.data %>>%
    dplyr::select(-matches('^true_')) %>>%
    dplyr::mutate(rm=rowMeans(dplyr::select(., -cancer))) %>>%
    ggplot(aes(rm, fill=cancer))+
    geom_histogram(position='dodge')+
    wtl::theme_wtl()

.data %>>%
    ggplot(aes(true_mu))+
    geom_histogram(aes(fill=cancer), bins=32, position='dodge')+
    wtl::theme_wtl()

.geno_matrices = .data %>>%
    dplyr::select(-matches('^true_')) %>>%
    tidyr::nest(-cancer) %>>%
    dplyr::mutate(data=purrr::map(data, as.matrix)) %>>% (?.)

.genotype_cancer = dplyr::filter(.geno_matrices, cancer)$data[[1]] %>>%
    (?as_tibble(.))

.genotype_cancer %>>% wtl::write_df('genotype1.tsv')

#########1#########2#########3#########4#########5#########6#########7#########
# interaction terms

genotypes2 = make_genotypes2(genotypes) %>>% (?.)
true_coefs2 = genotypes2 %>>%
    {setNames(runif(ncol(.)), names(.))} %>>%
    (. * (1.0 - true_epsilon) / sum(.)) %>>%
    (?.) %>>% (?sum(.))
stopifnot(sum(true_coefs2, true_epsilon) == 1.0)

.data2 = make_phenotypes(as.matrix(genotypes2), true_coefs2, sd_coef, true_threshold, true_epsilon) %>>%
    (?.) %>>% (? dplyr::count(., cancer))

.genotype_cancer2 = .data2 %>>%
    dplyr::filter(cancer) %>>%
    dplyr::select(-matches('^true_|^cancer$|:')) %>>% (?.)

.genotype_cancer2 %>>% wtl::write_df('genotype2.tsv')

for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(.genotype_cancer2, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    wtl::write_df(.df, sprintf('genotype2_%s.tsv', .term))
}

#########1#########2#########3#########4#########5#########6#########7#########
# Shuffle genotypes

working_dir = '~/Dropbox/working/tumor'
.infile = file.path(working_dir, 'genotype_pereira+er.tsv')
genotype_pereira = readr::read_tsv(.infile)
genotype_pereira %>>% dplyr::summarise_all(mean)
genotype_pereira %>>% as.matrix() %>>% table()
genotype_pereira %>>% dplyr::select(-Other) %>>% as.matrix() %>>% table()

genotype_pereira %>>%
    dplyr::mutate(rm=rowMeans(.)) %>>%
    ggplot(aes(rm))+geom_histogram()+wtl::theme_wtl()

shuffle_genotypes = function(.data) {
    .names = names(.data)
    purrr::by_row(.data, .labels=FALSE,
        ~{sample(.x) %>>% setNames(.names)}) %>>%
    tidyr::unnest()
}

.shuffled = genotype_pereira %>>% shuffle_genotypes() %>>% (?.)
.outfile = file.path(working_dir, sprintf('genotype_pereira+er_%s.tsv', wtl::now()))
wtl::write_df(.shuffled, .outfile)

.dimensions = ncol(genotype_pereira)
.combn = combn(.dimensions, 2)
for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(genotype_pereira, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    wtl::write_df(.df, sprintf('pereira+er_%s.tsv', .term))
}

for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(.shuffled, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    wtl::write_df(.df, sprintf('shuffled+er_%s.tsv', .term))
}

#########1#########2#########3#########4
# Run C++

system('./a.out genotype2.tsv')

.infiles = list.files('genotypes', pattern='.+\\.tsv', full.names=TRUE) %>>% (?.)
.outfiles = str_replace(.infiles, '^genotypes/', 'loglik/') %>>% (?.)
.commands = paste('./a.out -c0.3 -e0.1 -g91 -n1000 -o', .outfiles, .infiles) %>>% (?.)

library(doParallel)

system_par = function(.commands, .cores=parallel::detectCores(logical=FALSE)) {
    cluster = parallel::makeCluster(.cores)
    on.exit(parallel::stopCluster(cluster))
    doParallel::registerDoParallel(cluster)
    foreach::foreach(command=.commands, .combine=c) %dopar% {
        system(command, ignore.stdout=TRUE, ignore.stderr=TRUE)
    }
}
system_par(.commands)

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
