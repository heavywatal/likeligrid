library(pipeR)
library(stringr)
library(tidyverse)
library(jsonlite)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())

repo = '~/git/innanlab'
cachedir = '~/Dropbox/working/cancer'
raw_maf = file.path(cachedir, 'drivers.primary.PASS.extracted.tsv.gz') %>>% read_tsv() %>>%
    dplyr::select(-txid) %>>% (?.)

.ncases = raw_maf %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::count(cancer_type, sort=TRUE) %>>% (?.)

(.cancer_types = dplyr::filter(.ncases, nn > 400L)$cancer_type)

## TODO
# Assumption of impact?
# Count gene only once?
maf = raw_maf %>>%
    dplyr::filter(impact != 'LOW') %>>%
    dplyr::distinct(cancer_type, sample, symbol) %>>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>>%
    dplyr::arrange(cancer_type) %>>% (?.)

maf %>>% dplyr::filter(cancer_type %in% .cancer_types) %>>% dplyr::count(sample) %>>% (quantile(.$n, c(0.5, 0.95, 0.99, 1)))

(.hypermut_quantile = maf %>>% dplyr::count(sample) %>>% (quantile(.$n, c(0.5, 0.95, 0.99, 1))))
.hypermut_quantile['99%']
(.hypermut_df = tibble(p=names(.hypermut_quantile), x=.hypermut_quantile) %>>% dplyr::filter(p != '100%'))

.p = maf %>>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::ungroup() %>>%
    ggplot(aes(n))+
    geom_histogram(bins=40)+
    geom_vline(data=.hypermut_df, aes(xintercept=x, colour=p))+
    scale_colour_manual(values=c('50%'='black', '95%'='brown', '99%'='tomato'))+
    facet_wrap(~cancer_type)+
    scale_x_log10(breaks=wtl::breaks_log10, labels=wtl::labels_log10)+
    labs(x='# mutated sites', y='# patients')+
    theme(legend.position=c(0.6, 0.05), legend.justification=c(0, 0), legend.direction='horizontal')
.p
ggsave('hist_mutated_sites.png', .p, width=12, height=7)

#########1#########2#########3#########4#########5#########6#########7#########

pathdefs = read_json('pathway-defs.json', simplifyVector=TRUE)
lengths(pathdefs)
pathdefs = pathdefs[c('vogelstein', 'vogelstein3', 'pereira', 'HALLMARK')]

join_pathway = function(.genotypes, .pathdef) {
    .pd = map_df(.pathdef, ~tibble(symbol=.x), .id='pathway')
    dplyr::left_join(.genotypes, .pd, by='symbol') %>>%
    dplyr::mutate(pathway= dplyr::coalesce(pathway, 'other'))
}
join_pathway(head(maf), pathdefs[['HALLMARK']])

.joined = purrr::map_df(pathdefs, .id='definition', ~join_pathway(maf, .x))

.tidy = .joined %>>%
    dplyr::distinct(cancer_type, sample, symbol, definition, pathway) %>>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>>%
    dplyr::arrange(cancer_type, sample, symbol) %>>%
    (?.)

.counts = dplyr::full_join(by=c('cancer_type', 'sample'),
    maf %>>%
      dplyr::count(cancer_type, sample) %>>%
      dplyr::ungroup() %>>% dplyr::rename(genes=n),
    .tidy %>>%
      dplyr::count(cancer_type, sample, definition, pathway) %>>%
      dplyr::count(cancer_type, sample, definition) %>>%
      dplyr::ungroup() %>>% dplyr::rename(pathways=nn))

.plot = function(.data, .title='', .max_g=max(.data$genes), .max_p=max(.data$pathways)) {
    .scatter = .data %>>%
        ggplot(aes(genes, pathways))+
        geom_jitter(width=0, height=0.4, alpha=0.25)+
        geom_vline(data=.hypermut_df, aes(xintercept=x, colour=p))+
        scale_colour_manual(values=c('50%'='black', '95%'='brown', '99%'='tomato'))+
        coord_cartesian(xlim=c(0, .max_g), ylim=c(0, .max_p))+
        theme(legend.position='none')
    .bar_g = .data %>>%
        ggplot(aes(genes))+
        geom_histogram(bins=30)+
        labs(title=.title)+
        coord_cartesian(xlim=c(0, .max_g))+
        theme(axis.title.x=element_blank())
    .bar_p = .data %>>%
        ggplot(aes(pathways))+
        geom_bar()+
        coord_flip(xlim=c(0, .max_p))+
        theme(axis.title.y=element_blank())
    cowplot::plot_grid(.bar_g, NULL, .scatter, .bar_p, rel_widths=c(2, 1), rel_heights=c(1, 2))
}

plot_mutdist = function(.data, .definition='HALLMARK', .quantile='99%', .n=12) {
    .max_g = .hypermut_quantile[.quantile]
    .data %>>%
    dplyr::filter(definition == .definition, genes <= .max_g) %>>%
    tidyr::nest(-cancer_type, -definition) %>>%
    head(.n) %>>%
    dplyr::mutate(plt= purrr::map2(data, paste(definition, cancer_type), .plot,
        .max_g=.max_g, .max_p=length(pathdefs[[.definition]]))) %>>%
    (cowplot::plot_grid(plotlist=.$plt, ncol=4))
}
plot_mutdist(.counts, .n=4)

tidyr::crossing(def= names(pathdefs), p= c('95%', '99%', '100%')) %>>%
purrr::by_row(~{
    .outfile = sprintf('mutdist-%s-%03d.png', .$def, readr::parse_number(.$p))
    message(.outfile)
    .p = plot_mutdist(.counts, .$def, .$p)
    ggsave(.outfile, .p, width=16, height=10)
})

(hypermutants = .counts %>>% dplyr::filter(genes > .hypermut_quantile['99%']) %>>% (sample) %>>% unique())

.matrices = .tidy %>>%
    dplyr::filter(cancer_type %in% .cancer_types, !sample %in% hypermutants) %>>%
    tidyr::nest(-definition) %>>%
    dplyr::mutate(data= purrr::map(data, ~{
        .x %>>%
        dplyr::count(cancer_type, sample, pathway) %>>%
        dplyr::ungroup() %>>%
        tidyr::spread(pathway, n, fill=0L) %>>%
        dplyr::select(-sample) %>>%
        tidyr::nest(-cancer_type)
    })) %>>%
    tidyr::unnest() %>>% (?.)

.plot = function(.data, .name='') {
    .data %>>%
    (tibble(pathway=names(.), s=colSums(.))) %>>%
    ggplot(aes(pathway, s))+
    geom_col()+
    labs(title=.name)+
    wtl::theme_wtl(base_size=8)+
    theme(axis.ticks=element_blank(), panel.grid.major.x=element_blank())
}
.plot(.matrices$data[[1]])

.p = .matrices %>>%
    dplyr::arrange(definition, cancer_type) %>>%
    dplyr::mutate(plt= purrr::map2(data, paste(definition, cancer_type), .plot)) %>>%
    (plt) %>>% (cowplot::plot_grid(plotlist=., ncol=12))
.p
ggsave('pathway_freqs.png', .p, width=18, height=9)
