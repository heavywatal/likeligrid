library(pipeR)
library(stringr)
library(tidyverse)
library(jsonlite)
library(wtl)

repo = '~/git/innanlab'
cachedir = '~/Dropbox/working/cancer'
maf = file.path(cachedir, 'drivers.primary.PASS.extracted.tsv.gz') %>>% read_tsv() %>>%
    dplyr::filter(impact != 'LOW') %>>%
    dplyr::select(-txid) %>>% (?.)

pathdefs = read_json('pathway-defs.json', simplifyVector=TRUE)
lengths(pathdefs)
.defnames = c('vogelstein', 'vogelstein3', 'pereira', 'HALLMARK')

.ncases = maf %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::count(cancer_type, sort=TRUE) %>>% (?.)

(.cancer_types = dplyr::filter(.ncases, nn > 400L)$cancer_type)

.nested = maf %>>%
    tidyr::nest(-cancer_type) %>>%
    dplyr::filter(cancer_type %in% .cancer_types) %>>%
    (setNames(.$data, .$cancer_type)) %>>% (?.)

.maf = .nested[['BRCA']]
.pdef = pathdefs[['HALLMARK']]

make_matrix = function(.genotypes, .pathdef) {
    .pd = map_df(.pathdef, ~tibble(symbol=.x), .id='pathway')
    .genotypes %>>%
    dplyr::left_join(.pd, by='symbol') %>>%
    dplyr::count(sample, pathway) %>>%
    dplyr::ungroup() %>>%
    dplyr::mutate(pathway= dplyr::coalesce(pathway, 'other')) %>>%
    tidyr::spread(pathway, n, fill=0L) %>>%
    dplyr::select(-sample)
}
make_matrix(.maf, .pdef)

.matrices = purrr::map_df(pathdefs[.defnames], .id='definition', ~{
    .pd = .x
    purrr::map(.nested, ~{
        make_matrix(.x, .pd)
    }) %>>% (tibble(definition=names(.nested), data=.))
}) %>>% (?.)


.plot = function(.data) {
    .data %>>%
    (tibble(pathway=names(.), s=colSums(.))) %>>%
    ggplot(aes(pathway, s))+
    geom_col()+
    wtl::theme_wtl()
}
.plot(.matrices$data[[1]])

.matrices %>>%
    dplyr::mutate(plt= purrr::map(data, .plot))
