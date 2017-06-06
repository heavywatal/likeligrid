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

.plot = function(uniaxis, limits, gradient=tibble(), label='', ...) {
    .title = label
    .mle = dplyr::top_n(uniaxis, 1L, loglik) %>% dplyr::distinct(loglik, pathway, .keep_all=TRUE)
    .max_loglik = .mle$loglik[1]
    .p = ggplot(uniaxis, aes(value, loglik))+
        geom_line(colour='#66cc66')+
        geom_point(data=limits, alpha=0.5)
    .max_y = .max_loglik
    if (nrow(gradient) > 0L) {
        .offset = .max_loglik - min(gradient$loglik)
        .data_g = gradient %>%
          tidyr::gather(pathway, value, -loglik) %>%
          dplyr::distinct() %>%
          dplyr::mutate(loglik=loglik+.offset)
        .p = .p + geom_point(data=.data_g, colour='#33cc99', alpha=0.5)
        .max_y = max(.data_g$loglik)
    }
    .p+
    geom_vline(data=.mle, aes(xintercept=value), colour='salmon')+
    geom_hline(data=.chisq, aes(yintercept=.max_loglik - y, colour=alpha))+
    coord_cartesian(ylim=c(.max_y, .max_loglik - qchisq(0.95, 1)))+
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

# #######1#########2#########3#########4#########5#########6#########7#########
# gradient descent

.results_with_g = .results %>%
    dplyr::filter(str_detect(group, 'vogelstein$')) %>%
    dplyr::filter(str_detect(label, '-s4$')) %>% print() %>%
    dplyr::mutate(
        gradient=map(indir, ~{
            .f = list.files(.x, 'gradient-s5\\.tsv\\.gz', full.names=TRUE)
            message(.f)
            if (length(.f) > 0) {read_likeligrid(.f)} else {tibble()}
        })
    ) %>% print()

.out_g = .results_with_g %>% dplyr::mutate(plt= purrr::pmap(., .plot))
ggsave('gradient5.pdf', .out_g$plt, width=5, height=5)

# #######1#########2#########3#########4#########5#########6#########7#########
# epistasis

.results_e = .results %>%
    dplyr::filter(str_detect(group, 'vogelstein$')) %>%
    dplyr::filter(str_detect(label, '-s4$')) %>%
    # dplyr::filter(str_detect(label, 'BRCA')) %>%
    dplyr::mutate(
        epistasis=map(indir, ~{
            message(.x)
            list.files(.x, 'gradient-s4-e.*\\.gz$', full.names=TRUE) %>%
            map(read_likeligrid)
        })
    ) %>% print()

.rename_tail = function(.names) {
    c(head(.names, -1L), paste0('~', tail(.names, 1L)))
}
.rename_tail(letters)

.out_e = .results_e %>% dplyr::mutate(plts= purrr::pmap(., function(...) {
    .args = list(...)
    message(.args$label)
    map(.args$epistasis, ~{
        names(.x) = .rename_tail(names(.x))
        .plot(.args$uniaxis, .args$limits, .x, .args$label)
    })
}))
cowplot::plot_grid(plotlist=.out_e$plts[[1]])
.pl = .out_e$plts %>% map(~{cowplot::plot_grid(plotlist=.x)})
ggsave('epistasis-less.pdf', .pl, width=20, height=20)
