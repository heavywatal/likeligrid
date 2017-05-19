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
        dplyr::summarise(s= n(), p= prod(p)) %>>%
        dplyr::left_join(.dups, by=c('sample', 'pathway')) %>>%
        dplyr::summarise(p= prod(p, coef_excl, na.rm=TRUE) * permut(s), s= sum(s))
}
calc_numer(erpos, .freqs, .excl)

#########1#########2#########3#########4#########5#########6#########7#########

# .freqs = c(a=0.5, b=0.5)
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
    purrr::pmap_dbl(function(...) {
        .x = list(...)
        dups = table(purrr::flatten_chr(.x)) - 1L
        prod(x[names(dups)] ^ dups)
    })
}
calc_coefs_excl(.excl, 3L)

calc_denom = function(weights, excl, times=seq_len(6L)) {
    purrr::map_df(times, ~{
        .d = sum(calc_probs_null(weights, .x) * calc_coefs_excl(excl, .x))
        tibble::tibble(s=.x, denom= .d)
    })
}
calc_denom(.freqs, .excl, seq_len(3L))


calc_loglik = function(.data, excl=0.6, weights=.freqs, max_s=4L) {
    .excl = setNames(purrr::rep_along(.freqs, excl), names(weights))
    calc_numer(.data, weights, .excl) %>>%
    dplyr::filter(s <= max_s) %>>%
    dplyr::left_join(calc_denom(weights, .excl, seq_len(max_s)), by='s') %>>%
    dplyr::summarise(loglik= sum(log(p / denom))) %>>% (loglik)
}

.results = tibble::tibble(excl= seq(0.1, 1.0, 0.1)) %>>%
    dplyr::mutate(loglik= map_par(excl, ~calc_loglik(erpos, .x, max_s=6L)) %>>% flatten_dbl()) %>>%
    (?.)

.results


#########1#########2#########3#########4#########5#########6#########7#########
## Test

.genotypes = 'A,B
2,0
2,0
0,2
1,1
1,1
1,1
1,1
1,1
1,1
1,1
' %>>% read_csv() %>>% (?.)

.genotypes %>>% write_tsv('test_genotypes.tsv')

.n = nrow(.genotypes)
.excl = c(0.5, 0.4)
.sums = colSums(.genotypes)
.x = .sums / sum(.sums)
.c = sum(.excl * .x ^ 2) + 2 * prod(.x)
.dups = .genotypes %>>% mutate_all(function(x) {ifelse(x > 0L, x - 1, x)}) %>>% colSums()
.polynoms = 7 * log(2)

sum(.dups * log(.excl) + .sums * log(.x)) - .n * log(.c) + .polynoms
