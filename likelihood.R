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

.indir = '~/Dropbox/working/cancer'
.infile = file.path(.indir, 'pereira_tidy.tsv')
pereira_tidy = read_tsv(.infile) %>>% (?.)

erpos = pereira_tidy %>>%
    dplyr::filter(erStatus == 'POS', pathway != 'Other') %>>%
    dplyr::select(sample, pathway) %>>%
    (?.)

erpos %>>% dplyr::count(sample) %>>% group_by(n) %>>% summarise(samples=n()) %>>% dplyr::arrange(n)

permut = function(x) {
    factorial(sum(x)) / prod(factorial(x))
}

calc_numer = function(.data, weights, excl) {
    .dups = .data %>>%
        dplyr::filter(duplicated(.)) %>>%
        dplyr::mutate(coef_excl= excl[pathway])
    .data %>>%
        dplyr::mutate(p= weights[pathway]) %>>%
        dplyr::group_by(sample, pathway) %>>%
        dplyr::summarise(s= n(), p= prod(p)) %>>% (?.) %>>%
        dplyr::left_join(.dups, by=c('sample', 'pathway')) %>>%
        dplyr::summarise(p= prod(p, coef_excl, na.rm=TRUE) * permut(s), s= sum(s))
}
calc_numer(erpos, .freqs, .excl)

#########1#########2#########3#########4#########5#########6#########7#########

.freqs = erpos %>>% dplyr::count(pathway) %>>% (setNames(.$n / sum(.$n), .$pathway))
.excl = setNames(purrr::rep_along(.freqs, 0.6), names(.freqs))

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

calc_denom = function(weights, excl, times=seq_len(6L)) {
    purrr::map_df(times, ~{
        .d = sum(calc_probs_null(weights, .x) * calc_coefs_excl(excl, .x))
        tibble::tibble(s=.x, denom= .d)
    })
}
calc_denom(.freqs, .excl, seq_len(3L))
