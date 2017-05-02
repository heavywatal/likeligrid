library(pipeR)
library(tidyverse)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())
#########1#########2#########3#########4#########5#########6#########7#########

read_likeligrid = function(file) {
    read_tsv(file, comment='#', col_types=cols(.default='d'))
}

.extract_mle = function(.data) {
    dplyr::top_n(.data, 1, wt=loglik) %>>% dplyr::distinct() %>>% head(1)
}

.transform = function(.data, .best) {
    tibble(pathway= names(.data)[-1]) %>>%
    purrr::by_row(~{
        .data[.data[.x$pathway] != .best[[.x$pathway]],] %>>%
        dplyr::bind_rows(.best) %>>%
        dplyr::arrange_(.x$pathway) %>>%
        dplyr::transmute_(value= .x$pathway, ~loglik)
    }) %>>%
    tidyr::unnest()
}

.chisq = tibble(alpha=c(0.9, 0.95, 0.99), y=qchisq(alpha, 1) * 0.5)

.plot_tidy = function(.tidy, .best, .title='') {
    ggplot(.tidy, aes(value, loglik))+
    geom_line(size=rel(1.5))+
    geom_vline(aes(xintercept=x), data=.best %>>% tidyr::gather(pathway, x, -loglik), colour='salmon')+
    geom_hline(aes(yintercept=.best$loglik - y, colour=alpha), data=.chisq)+
    scale_colour_gradient(breaks=.chisq$alpha, guide='legend')+
    scale_x_continuous(breaks=c(0, 1, 2), limits=c(0, 2), expand=c(0, 0))+
    facet_wrap(~pathway)+
    labs(title=.title)+
    wtl::theme_wtl()
}

.plot = function(.data, .title='') {
    .best = .data %>>% .extract_mle()
    .data %>>%
        dplyr::filter(loglik > .best$loglik - qchisq(0.99, 1)) %>>%
        .transform(.best) %>>%
        .plot_tidy(.best, .title)
}

#########1#########2#########3#########4#########5#########6#########7#########

.outdir = '~/working/likeligrid/analysis'
# .datadir = '~/working/likeligrid/output-tiny-0501'
.datadir = '~/working/likeligrid/output-major-0501'

.uniaxis = tibble(indir= list.dirs(.datadir, recursive=FALSE)) %>>%
    dplyr::transmute(
        label= basename(indir),
        group= str_replace(label, '-s\\d$', ''),
        data= purrr::map(indir, ~{
            list.files(.x, 'uniaxis.+\\.gz', full.names=TRUE) %>>%
            purrr::map_df(read_likeligrid)
        })
    ) %>>% (?.)

.plot(.uniaxis$data[[2]], 'test')

.uniaxis %>>%
    dplyr::filter(purrr::map_int(data, nrow) > 0) %>>%
    dplyr::group_by(group) %>>%
    purrr::by_slice(~{
        purrr::by_row(.x, ~.plot(.x$data[[1]], .x$label)) %>>%
        (cowplot::plot_grid(plotlist=.$.out, nrow=2, ncol=2))
    }) %>>%
    # (ggsave('uniaxis-tiny-0501.pdf', .$.out, width=12, height=12))
    (ggsave('uniaxis-major-0501.pdf', .$.out, width=12, height=12))

#########1#########2#########3#########4#########5#########6#########7#########
# Plot including limit-*.tsv.gz

.read_limit_max = function(.infile) {
    .pathway = str_extract(.infile, '(?<=limit-)[^.]+(?=_[LU]\\.tsv)')
    read_likeligrid(.infile) %>>%
    dplyr::transmute_(.dots=c(pathway=~.pathway, value=.pathway, 'loglik')) %>>%
    dplyr::group_by(value) %>>%
    dplyr::top_n(1, wt=loglik) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange(value)
}

.gather_uniaxis = function(.data) {
    .mle = .extract_mle(.data)
    tibble(pathway= names(.data)[-1]) %>>%
    purrr::by_row(~{
        .data[.data[.x$pathway] != .mle[[.x$pathway]],] %>>%
        dplyr::bind_rows(.mle) %>>%
        dplyr::arrange_(.x$pathway) %>>%
        dplyr::transmute_(value= .x$pathway, ~loglik)
    }) %>>%
    tidyr::unnest()
}

.indirs = list.dirs(.datadir, recursive=FALSE)
.infiles = list.files(.indirs[15], 'limit.*\\.tsv\\.gz$', full.names=TRUE)
# .read_limit_max(.infiles[1])
.limits = .infiles %>>% map_df(.read_limit_max) %>>% (?.)
.uniaxis = file.path(.indirs[15], 'uniaxis.tsv.gz') %>>% read_likeligrid() %>>% (?.)

.plot = function(.limits, .uniaxis, .title='') {
    .mle = .extract_mle(.uniaxis)
    ggplot(.limits, aes(value, loglik))+
    geom_line(data=.gather_uniaxis(.uniaxis), colour='#66cc66')+
    geom_point(alpha=0.5)+
    geom_vline(data=.mle %>>% tidyr::gather(pathway, x, -loglik), aes(xintercept=x), colour='salmon')+
    geom_hline(data=.chisq, aes(yintercept=.mle$loglik - y, colour=alpha))+
    coord_cartesian(ylim=c(.mle$loglik, .mle$loglik- qchisq(0.95, 1)))+
    scale_x_continuous(breaks=c(0, 1, 2), limits=c(0, 2), expand=c(0, 0))+
    labs(title=.title)+
    facet_wrap(~pathway)+
    wtl::theme_wtl()+
    theme(legend.position='none', panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(colour='#aaaaaa'),
      axis.ticks=element_blank(), axis.title.x=element_blank())
}
.plot(.limits, .uniaxis)

plot_limits = function(.indir) {
    message(.indir)
    .infiles = list.files(.indir, 'limit.*\\.tsv\\.gz$', full.names=TRUE)# %>>% (?.)
    .limits = .infiles %>>% map_df(.read_limit_max)# %>>% (?.)
    .uniaxis = file.path(.indir, 'uniaxis.tsv.gz') %>>% read_likeligrid()# %>>% (?.)
    .plot(.limits, .uniaxis, basename(.indir))
}
plot_limits(.indirs[15])

.out = tibble(
    indir= list.dirs(.datadir, recursive=FALSE),
    group= str_replace(indir, '-s\\d$', '') %>>% basename()) %>>%
    nest(-group) %>>%
    # dplyr::filter(str_detect(group, '3$')) %>>% head(4) %>>%
    dplyr::mutate(plt= wtl::map_par(data, ~{
        .plts = purrr::map(.x$indir, plot_limits)
        cowplot::plot_grid(plotlist=.plts, nrow=2, ncol=2)
    }))

# .out$plt[[2]]
ggsave('confidence-simple.pdf', .out$plt, width=12, height=10)


#########1#########2#########3#########4#########5#########6#########7#########
# limits in detail

.indirs
.infiles = list.files(.indirs[4], 'limit.*\\.tsv\\.gz$', full.names=TRUE)

.uniaxis = file.path(.indirs[4], 'uniaxis.tsv.gz') %>>% read_likeligrid() %>>% (?.)
.mle = .extract_mle(.uniaxis)

.read_limit_max(.infiles[1])

.read_limit_max = function(.infile) {
    .pathway = str_extract(.infile, '(?<=limit-)[^.]+(?=_[LU]\\.tsv)')
    read_likeligrid(.infile) %>>%
    dplyr::group_by_(.pathway) %>>%
    dplyr::top_n(1, wt=loglik) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange_(.pathway)
}
.read_limit_max(.infiles[1])


.data = read_likeligrid(.infiles[1])

.data %>>%
    dplyr::group_by(Fate) %>>%
    dplyr::top_n(1, wt=loglik)

.mle
.data %>>% dplyr::summarise_all(max)
.data %>>% dplyr::summarise_all(min)
.data %>>% dplyr::arrange(desc(loglik))

dplyr::top_n(1, wt = loglik)
