library(tidyverse)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())
#########1#########2#########3#########4#########5#########6#########7#########

read_likeligrid = function(file) {
    read_tsv(file, comment='#', col_types=cols(.default='d'))
}

read_limit_max = function(.infile) {
    .pathway = str_extract(.infile, '(?<=limit-)[^.]+(?=_[LU]\\.tsv)')
    read_likeligrid(.infile) %>%
    dplyr::transmute_(.dots=c(pathway=~.pathway, value=.pathway, 'loglik')) %>%
    dplyr::group_by(value) %>%
    dplyr::top_n(1, wt=loglik) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(value)
}

gather_uniaxis = function(.data) {
    .data %>%
    dplyr::select_if(function(.x) {dplyr::n_distinct(.x) > 1 || .x[1] < 0}) %>%
    tidyr::gather(pathway, value, -loglik)
}

.chisq = tibble(alpha=c(0.9, 0.95, 0.99), y=qchisq(alpha, 1) * 0.5)

.plot = function(uniaxis, limits, label='') {
    .title = label
    .mle = dplyr::top_n(uniaxis, 1L, loglik) %>% dplyr::distinct(loglik, pathway, .keep_all=TRUE)
    .max_loglik = .mle$loglik[1]
    ggplot(limits, aes(value, loglik))+
    geom_line(data=uniaxis, colour='#66cc66')+
    geom_point(alpha=0.5)+
    geom_vline(data=.mle, aes(xintercept=value), colour='salmon')+
    geom_hline(data=.chisq, aes(yintercept=.max_loglik - y, colour=alpha))+
    coord_cartesian(ylim=c(.max_loglik, .max_loglik - qchisq(0.95, 1)))+
    scale_x_continuous(breaks=c(0, 1, 2), limits=c(0, 2), expand=c(0, 0))+
    labs(title=.title)+
    facet_wrap(~pathway)+
    wtl::theme_wtl()+
    theme(legend.position='none', panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(colour='#aaaaaa'),
      axis.ticks=element_blank(), axis.title.x=element_blank())
}

.datadir = '~/working/likeligrid/output'
.result_dirs = tibble(indir= list.dirs(.datadir, recursive=FALSE)) %>%
    dplyr::mutate(
        label= basename(indir),
        group= str_replace(label, '-s\\d$', '')) %>%
    dplyr::filter(str_detect(group, 'vogelstein$')) %>% print()

.results = .result_dirs %>%
    dplyr::filter(map_lgl(indir, ~{
        list.files(.x, 'limit.+\\.gz$') %>% length() %>% {. > 0L}
    })) %>%
    dplyr::mutate(
        uniaxis= purrr::map(indir, ~{
            list.files(.x, 'uniaxis.+\\.gz$', full.names=TRUE) %>%
            purrr::map_df(~{read_likeligrid(.x) %>% gather_uniaxis()})
        }),
        limits= purrr::map(indir, ~{
            list.files(.x, 'limit.+\\.gz$', full.names=TRUE) %>%
            purrr::map_df(read_limit_max)
        })
    ) %>% print()

.out = .results %>%
    dplyr::select(group, label, uniaxis, limits) %>%
    tidyr::nest(-group) %>%
    dplyr::mutate(plt= purrr::map(data, ~{
        .plts = .x %>%
            dplyr::filter(purrr::map_lgl(limits, ~{nrow(.x) > 0})) %>%
            {message(.$label); .} %>%
            purrr::pmap(.plot)
        cowplot::plot_grid(plotlist=.plts, nrow=2, ncol=2)
    }))

# .out$plt[[2]]
ggsave('confidence.pdf', .out$plt, width=9, height=9)
