library(pipeR)
library(tidyverse)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())
#########1#########2#########3#########4#########5#########6#########7#########

.filter = function(.data, .best, .column) {
    .data[.data[.column] != .best[[.column]],]
}

.extract_best = function(.data) {
    .data[.data %>>% duplicated(),][1,]
}

.transform = function(.data, .best) {
    tibble(pathway= names(.data)[-1]) %>>%
    purrr::by_row(~{
        .filter(.data, .best, .x$pathway) %>>%
        dplyr::bind_rows(.best) %>>%
        dplyr::arrange_(.x$pathway) %>>%
        dplyr::transmute_(value= .x$pathway, ~loglik)
    }) %>>%
    tidyr::unnest()
}

.chisq = tibble(alpha=c(0.9, 0.95, 0.99), y=qchisq(alpha, 1) * 0.5)

.plot = function(.data, .best, .title='') {
    ggplot(.data, aes(value, loglik))+
    geom_line(size=rel(1.5))+
    geom_vline(aes(xintercept=x), data=.best %>>% tidyr::gather(pathway, x, -loglik), colour='salmon')+
    geom_hline(aes(yintercept=.best$loglik - y, colour=alpha), data=.chisq)+
    scale_colour_gradient(breaks=.chisq$alpha, guide='legend')+
    scale_x_continuous(breaks=c(0, 1, 2), limits=c(0, 2), expand=c(0, 0))+
    facet_wrap(~pathway)+
    labs(title=.title)+
    wtl::theme_wtl()
}

plot_uniaxis = function(.infiles) {
    tibble(infile=.infiles, label=basename(dirname(infile))) %>>%
    purrr::by_row(~{
        .data = read_tsv(.x$infile, comment='#')
        .best = .extract_best(.data)
        .data = .data %>>% dplyr::filter(loglik > .best$loglik - qchisq(0.99, 1))
        .transform(.data, .best) %>>%
        .plot(.best, .x$label)
    }) %>>%
    (cowplot::plot_grid(plotlist=.$.out, nrow=2, ncol=2))
}
plot_uniaxis(.infiles)

#########1#########2#########3#########4#########5#########6#########7#########

.outdir = '~/working/likeligrid/analysis'
.datadir = '~/working/likeligrid/output'

.nested = tibble(indir= list.dirs(.datadir, recursive=FALSE)) %>>%
    dplyr::transmute(
        infile= file.path(indir, 'uniaxis.tsv.gz'),
        group= str_replace(indir, '-s\\d$', '') %>>% basename()
    ) %>>% tidyr::nest(-group) %>>%
    dplyr::mutate(data= purrr::map(data, flatten_chr)) %>>% (?.)

.plts = .nested %>>%
    dplyr::mutate(plts= wtl::map_par(data, plot_uniaxis)) %>>%
    (plts)

grid.draw.list = function(x, ...) {
    purrr::walk(x, grid::grid.draw, ...)
}
# grid::grid.draw(.plts)

ggsave('uniaxis.pdf', .plts, width=12, height=12)
