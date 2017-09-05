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
    dplyr::select(loglik, one_of(.pathway)) %>%
    setNames(c('loglik', 'value')) %>%
    dplyr::mutate(pathway=.pathway) %>%
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
    if (gradient %>% inherits('list')) {
        return(purrr::map(gradient, ~.plot(uniaxis, limits, .x, label, ...)))
    }
    .mle = dplyr::top_n(uniaxis, 1L, loglik) %>% dplyr::distinct(loglik, pathway, .keep_all=TRUE)
    .max_loglik = .mle$loglik[1]
    .p = ggplot(uniaxis, aes(value, loglik))+
        geom_line(colour='#66cc66')+
        geom_point(data=limits, alpha=0.5)
    .max_y = .max_loglik
    if (nrow(gradient) > 0L) {
        if (nrow(gradient) > 100L) {
            gradient = gradient %>%
                dplyr::arrange(desc(loglik)) %>%
                {dplyr::bind_rows(head(., 3L), dplyr::sample_n(., 80L), tail(., 3L))} %>%
                dplyr::distinct()
        }
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
    labs(title=label)+
    facet_wrap(~pathway)+
    wtl::theme_wtl()+
    theme(legend.position='none', panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(colour='#aaaaaa'),
      axis.ticks=element_blank(), axis.title.x=element_blank())
}

.datadir = '~/working/likeligrid/output'

.metadata = tibble(indir= list.dirs(.datadir, recursive=FALSE)) %>%
    dplyr::mutate(label= basename(indir)) %>%
    tidyr::separate(label, c('TCGA', 'type', 'definition', 's', 'epistasis'), '-', fill='right', remove=FALSE) %>%
    dplyr::select(-TCGA) %>%
    dplyr::filter(!is.na(type)) %>%
    dplyr::mutate(s=str_replace(s, '^s', '') %>% parse_integer()) %>%
    print()

.results = .metadata %>%
    dplyr::filter(is.na(epistasis)) %>%
    dplyr::filter(s == 4L) %>%
    dplyr::mutate(
        uniaxis= purrr::map(indir, ~{
            message(.x)
            list.files(.x, 'uniaxis.+\\.gz$', full.names=TRUE) %>%
            purrr::map_dfr(~{read_likeligrid(.x) %>% gather_uniaxis()})
        }),
        limits= purrr::map(indir, ~{
            message(.x)
            list.files(.x, 'limit.+\\.gz$', full.names=TRUE) %>%
            purrr::map_dfr(read_limit_max)
        })
    ) %>% print()

.results %>%
    tidyr::nest(-type) %>%
    dplyr::mutate(plt= purrr::map(data, ~{
        .plts = .x %>%
            dplyr::filter(purrr::map_lgl(limits, ~{nrow(.x) > 0})) %>%
            {message(.$label); .} %>%
            purrr::pmap(.plot)
        cowplot::plot_grid(plotlist=.plts, nrow=2, ncol=2)
    }))

# .out$plt[[2]]
ggsave('confidence.pdf', .out$plt, width=9, height=9)

.plt = .results %>% purrr::pmap(.plot)
.p = cowplot::plot_grid(plotlist=.plt, nrow=3, ncol=4)
ggsave('confidence-s4.pdf', .p, width=12, height=9)

# #######1#########2#########3#########4#########5#########6#########7#########
# gradient descent

.results_with_g = .results %>%
    dplyr::filter(str_detect(definition, 'vogelstein$')) %>%
    dplyr::filter(str_detect(label, '-s4$')) %>% print() %>%
    dplyr::mutate(
        gradient=purrr::map(indir, ~{
            .f = list.files(.x, 'gradient-s5\\.tsv\\.gz', full.names=TRUE)
        })
    ) %>% print()

.plt = .results_with_g %>% purrr::pmap(.plot)
.p = cowplot::plot_grid(plotlist=.plt, nrow=3, ncol=4)
ggsave('gradient5.pdf', .p, width=12, height=9)
# ggsave('gradient6.pdf', .p, width=12, height=9)

# #######1#########2#########3#########4#########5#########6#########7#########
# epistasis

# Put interaction term at last of facet_wrap
.rename_tail = function(.names) {
    c(head(.names, -1L), paste0('~', tail(.names, 1L)))
}
.rename_tail(letters)

.results_e = .results %>%
    dplyr::filter(str_detect(definition, 'vogelstein$')) %>%
    dplyr::filter(str_detect(label, '-s4$')) %>%
    # dplyr::filter(str_detect(label, 'BRCA')) %>%
    dplyr::mutate(
        gradient=purrr::map(indir, ~{
            message(.x)
            list.files(.x, 'gradient-s5-e.*\\.gz$', full.names=TRUE) %>%
            purrr::map(~{
                read_likeligrid(.x) %>%
                setNames(names(.) %>% .rename_tail())
            })
        })
    ) %>% print()

.results_e$gradient[[1]] %>% map_int(nrow)

.plts_epistasis = .results_e %>% purrr::pmap(.plot)
.pl = .plts_epistasis %>% purrr::map(~{cowplot::plot_grid(plotlist=.x)})
# ggsave('epistasis-s4.pdf', .pl, width=20, height=20)
ggsave('epistasis-s5.pdf', .pl, width=20, height=20)
