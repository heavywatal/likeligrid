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
# Shuffle genotypes

working_dir = '~/Dropbox/working/cancer/genotypes'
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

shuffle_genotype_file = function(.infile) {
    .data = readr::read_tsv(.infile)
    .outfile = str_replace(.infile, '\\bgenotype_', 'shuffled_')
    message(.outfile)
    shuffle_genotypes(.data) %>>% write_tsv(.outfile, na='')
}

working_dir = '~/Dropbox/working/cancer/genotypes'
.genotype_files = list.files(working_dir, pattern='\\bgenotype_pereira', full.names=TRUE)
shuffle_genotype_file(.genotype_files[1])

.genotype_files %>>% purrr::walk(shuffle_genotype_file)

.shuffled = genotype_pereira %>>% shuffle_genotypes() %>>% (?.)
.outfile = file.path(working_dir, sprintf('genotype_shuffled+er_%s.tsv', wtl::now()))
write_tsv(.shuffled, .outfile, na='')

.dimensions = ncol(genotype_pereira)
.combn = combn(.dimensions, 2)
for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(genotype_pereira, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    write_tsv(.df, sprintf('pereira+er_%s.tsv', .term), na='')
}

for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(.shuffled, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    write_tsv(.df, sprintf('shuffled+er_%s.tsv', .term), na='')
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

.genotype_cancer %>>% write_tsv('genotype1.tsv', na='')

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

.genotype_cancer2 %>>% write_tsv('genotype2.tsv', na='')

for (i in seq_len(ncol(.combn))) {
    .df = add_interaction_term(.genotype_cancer2, .combn[,i])
    .term = tail(names(.df), 1)
    message(.term)
    write_tsv(.df, sprintf('genotype2_%s.tsv', .term), na='')
}

#########1#########2#########3#########4
# Run C++

system('./a.out genotype2.tsv')

.infiles = list.files('genotypes', pattern='.+\\.tsv', full.names=TRUE) %>>% (?.)
.outfiles = str_replace(.infiles, '^genotypes/', 'loglik/') %>>% (?.)
.commands = paste('./a.out -c0.3 -e0.1 -g21 -n1000 -o', .outfiles, .infiles) %>>% (?.)

wtl::map_par(.commands, system)
