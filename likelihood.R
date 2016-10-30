.indir = '~/Dropbox/working/cancer'
.infile = file.path(.indir, 'genotypes/genotype_pereira+er.tsv')
genotype_pereira = read_tsv(.infile) %>>% (?.)

.names = names(genotype_pereira)
.freqs = genotype_pereira %>>% colMeans() %>>% (. / sum(.))

.nested = genotype_pereira %>>%
    dplyr::group_by(s=rowSums(.)) %>>%
    tidyr::nest() %>>%
    dplyr::arrange(s) %>>%
    (?.)

.nested$data[[3]]

foreach::foreach(row=iterators::iter(.nested$data[[3]] %>>% head(), by='row')) %do% {
    .v = unlist(row)
    .freqs ^ .v
}

#########1#########2#########3#########4#########5#########6#########7#########

.infile = file.path(.indir, 'pereira_tidy.tsv')
pereira_tidy = read_tsv(.infile) %>>% (?.)

erpos = pereira_tidy %>>%
    dplyr::filter(erStatus == 'POS', pathway != 'Other') %>>%
    dplyr::select(sample, pathway) %>>%
    (?.)

permut = function(x) {
    factorial(sum(x)) / prod(factorial(x))
}

calc_loglik = function(.data, exclusivity=0.6) {.data %>>%
    dplyr::mutate(p= .freqs[pathway]) %>>%
    dplyr::group_by(sample) %>>%
    purrr::by_slice(~{
        s_paths = table(.$pathway)
        dups = sum(s_paths - 1L)
        log(prod(.$p) * permut(s_paths) * exclusivity^dups)
    }, .to= 'lnp') %>>%
    tidyr::unnest()
}
erpos %>>% calc_loglik()

seq(0.1, 1.0, 0.1) %>>%
    purrr::map_dbl(~{calc_loglik(erpos, .x)$lnp %>>% sum()})

#########1#########2#########3#########4#########5#########6#########7#########

.names = names(genotype_pereira)
.freqs = genotype_pereira %>>% colMeans() %>>% (. / sum(.))
.excl = setNames(purrr::rep_along(.names, 0.6), .names)

calc_probs_null = function(x, times) {
    crossing_rep(x, times) %>>%
    dplyr::transmute_(.out= paste(names(.), collapse='*')) %>>%
    `[[`('.out')
}
calc_probs_null(.freqs, 3L)

calc_coefs_excl = function(x, times) {
    crossing_rep(names(x), times) %>>%
    purrr::by_row(~{
        dups = table(purrr::flatten_chr(.x)) - 1L
        prod(x[names(dups)] ^ dups)
    }, .labels=FALSE) %>>%
    (flatten_dbl(.$.out))
}
calc_coefs_excl(.excl, 3L)

calc_denom = function(weights, excl, times) {
    sum(calc_probs_null(weights, times) * calc_coefs_excl(excl, times))
}
calc_denom(.freqs, .excl, 3L)
