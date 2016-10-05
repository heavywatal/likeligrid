library(pipeR)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(wtl)
#########1#########2#########3#########4#########5#########6#########7#########

setwd('~/working/cancer')
.drivers = read_tsv('driver_genes.tsv') %>>% (?.)
top_drivers = .drivers %>>%
    dplyr::filter(refs >= 4) %>>%
    dplyr::rename(pathway=vogel_pathway) %>>%
    tidyr::replace_na(list(role='TSG', pathway='NA')) %>>%
    dplyr::select(gene, role, pathway) %>>%
    (?.)
init_genes = setNames(rep(1, nrow(top_drivers)), top_drivers$gene) %>>% (?.)

.pathway_genes = .drivers %>>%
    tidyr::replace_na(list(vogel_pathway='NA')) %>>%
    dplyr::mutate(pathway= str_split(vogel_pathway, ':')) %>>%
    tidyr::unnest() %>>%
    dplyr::count(pathway) %>>%
    (?.)
init_pathways = setNames(rep(1, nrow(.pathway_genes)), .pathway_genes$pathway) %>>% (?.)

maf_annot = read_tsv('gdc_topdrivers.maf.clean.tsv.gz') %>>%
    dplyr::left_join(top_drivers, by='gene') %>>%
    (?.)

maf_annot %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::count(cancer_type) %>>%
    dplyr::arrange(desc(nn))

maf_annot %>>%
    dplyr::mutate(width = end - start) %>>%
    ggplot(aes(width))+geom_histogram(center=0, binwidth=1)

#########1#########2#########3#########4#########5#########6#########7#########

calc_mudens = function(.data) {.data %>>%
    tidyr::nest(-role, -pathway, -gene, -mutype) %>>%
    dplyr::select(-pathway) %>>%
    purrr::by_row(~{
        table(.$data[[1]]$start, dnn='pos') %>>%
        tibble::as_tibble() %>>%
        dplyr::mutate(pos=as.integer(pos),
          mut_density= n / sum(n), n=NULL)
    }) %>>%
    dplyr::select(-data) %>>%
    tidyr::unnest()
}

mut_density = maf_annot %>>%
    (? dplyr::count(., mutype)) %>>%
    dplyr::mutate(mutype= str_replace(mutype, 'splice', 'truncating')) %>>%
    dplyr::filter(mutype %in% c('missense', 'truncating')) %>>%
    calc_mudens() %>>% (?.)

mut_density %>>%
    dplyr::filter(mut_density > 0.1) %>>%
    dplyr::filter(!(role=='oncogene' & mutype=='truncating')) %>>%
    dplyr::arrange(role, gene, pos) %>>%
    (?.) %>>% wtl::write_df('recurrent_mutations.tsv.gz')

add_mudens_weight = function(.data) {
    .oncomut_density = mut_density %>>%
        dplyr::filter(role=='oncogene') %>>%
        dplyr::select(-role)
    left_join(.data, .oncomut_density, by=c('gene', 'pos', 'mutype')) %>>%
    dplyr::mutate(mudens_weight= ifelse(role=='oncogene', mut_density, 1)) %>>%
    tidyr::replace_na(list(mudens_weight=0)) %>>%
    dplyr::select(-mut_density)
}

#' @param tsg_missense numeric value relative to truncating mutation
summarise_mut_score_by_gene = function(.data, tsg_missense) {
    dplyr::mutate(.data, tsg_weight=
      ifelse(role!='oncogene' & mutype=='missense', tsg_missense, 1.0)) %>>%
    dplyr::group_by(sample, pathway, gene, role) %>>%
    dplyr::summarise(mut_score= sum(mudens_weight * tsg_weight)) %>>%
    dplyr::ungroup() %>>%
    dplyr::mutate(mut_score= pmin(1.0, mut_score))
}

maf_weighted = maf_annot %>>%
    add_mudens_weight() %>>%
    summarise_mut_score_by_gene(0.5) %>>%
    dplyr::arrange(sample, pathway, desc(mut_score)) %>>% (?.)

#########1#########2#########3#########4#########5#########6#########7#########

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
# Exhaustive loops

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
