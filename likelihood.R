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

exclusivity = 0.6

permut = function(x) {
    factorial(length(x)) / prod(factorial(table(x)))
}

erpos %>>%
    dplyr::mutate(p= .freqs[pathway]) %>>%
    dplyr::group_by(sample) %>>%
#    head(6) %>>%
    purrr::by_slice(~{
        log(prod(.$p) * permut(.$pathway))
    }, .to= 'lnp') %>>%
    tidyr::unnest()
