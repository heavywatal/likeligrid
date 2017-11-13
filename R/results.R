library(tidyverse)
library(wtl)
loadNamespace('cowplot')
#########1#########2#########3#########4#########5#########6#########7#########
source('~/git/likeligrid/R/common.R')

.datadir = '~/working/likeligrid/output'
.datadir = '~/working/likeligrid/output-pleiotropy'

.metadata = read_metadata(.datadir) %>% print()

# #######1#########2#########3#########4#########5#########6#########7#########

.results = .metadata %>%
    dplyr::filter(is.na(epistasis)) %>%
    dplyr::filter(s == 4L) %>%
    head(3) %>%
    dplyr::mutate(
        uniaxis= purrr::map(indir, read_all_uniaxis),
        limits= purrr::map(indir, read_all_limit)
    ) %>% print()

.plt = .results %>%
    # dplyr::filter(purrr::map_lgl(limits, ~{nrow(.x) > 0L})) %>%
    purrr::pmap(plot_landscape)
.p = cowplot::plot_grid(plotlist=.plt, nrow=3, ncol=4)
ggsave('confidence-s4.pdf', .p, width=12, height=9)

# #######1#########2#########3#########4#########5#########6#########7#########
# gradient descent

.results_with_g = .results %>%
    dplyr::filter(str_detect(definition, 'vogelstein$')) %>%
    dplyr::filter(s == 4L) %>%
    dplyr::mutate(
        gradient=purrr::map(indir, ~{
            list.files(.x, 'gradient-s5\\.tsv\\.gz', full.names=TRUE) %>%
            read_likeligrid()
        })
    ) %>% print()

.plt = .results_with_g %>% purrr::pmap(plot_landscape)
.p = cowplot::plot_grid(plotlist=.plt, nrow=3, ncol=4)
ggsave('gradient5.pdf', .p, width=12, height=9)
# ggsave('gradient6.pdf', .p, width=12, height=9)

# #######1#########2#########3#########4#########5#########6#########7#########
# epistasis

.results_e = .results %>%
    dplyr::filter(str_detect(definition, 'vogelstein$')) %>%
    dplyr::filter(s == 4L) %>%
    # dplyr::filter(str_detect(label, 'BRCA')) %>%
    dplyr::mutate(
        gradient=purrr::map(indir, ~{
            message(.x)
            list.files(.x, 'gradient-s5-e.*\\.gz$', full.names=TRUE) %>%
            purrr::map(read_likeligrid)
        })
    ) %>% print()

.results_e$gradient[[1]] %>% map_int(nrow)

.plts_epistasis = .results_e %>% purrr::pmap(plot_landscape)
.pl = .plts_epistasis %>% purrr::map(~{cowplot::plot_grid(plotlist=.x)})
# ggsave('epistasis-s4.pdf', .pl, width=20, height=20)
ggsave('epistasis-s5.pdf', .pl, width=20, height=20)

# #######1#########2#########3#########4#########5#########6#########7#########
# summarize likelihood elevation in gradient mode

.results_ge = .metadata %>%
    dplyr::filter(s == 4L) %>%
    dplyr::filter(is.na(epistasis)) %>%
    # head(3L) %>%
    dplyr::mutate(
      gradfiles = purrr::map(indir, ~list.files(.x, 'gradient-s5-e.*\\.gz$'))
    ) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      epistasis= str_extract(gradfiles, '(?<=-e)(\\dx\\d)'),
      grad_epi= purrr::map(file.path(indir, gradfiles), ~{
        message(.x)
        .df = read_likeligrid(.x) #%>% dplyr::arrange(loglik)
        .epistasis_col = .df %>% dplyr::select(ncol(.))
        .start = .df %>% dplyr::filter(.epistasis_col == 1.0) %>% top_n(1L, loglik)
        .end = .df %>% top_n(1L, loglik)
        dplyr::bind_rows(.start, .end) %>%
        tidyr::gather(pathway, value, -loglik)
      }),
      gradfiles= NULL
    ) %>%
    tidyr::unnest() %>%
    dplyr::mutate(pathway= str_replace(pathway, '\\w+:\\w+', 'epistasis')) %>%
    print()

.nested = .results_ge %>%
    dplyr::select(type, epistasis, loglik, pathway, value) %>%
    tidyr::nest(-type) %>%
    print()

.plot = function(.tbl, .title='') {
    .start_loglik = min(.tbl$loglik)
    ggplot(.tbl, aes(value, loglik))+
    geom_hline(data=.const.chisq, aes(yintercept=.start_loglik + y, colour=alpha))+
    geom_path(size=1.3)+
    facet_grid(epistasis ~ pathway)+
    coord_cartesian(xlim=c(0, 2))+
    scale_x_continuous(breaks=c(0, 1, 2))+
    labs(title=.title)+
    theme_wtl()+
    theme(legend.position='none')
}
.plot(.nested$data[[1]], 'test')

.plt = .nested %>%
    dplyr::mutate(plt= purrr::map2(data, type, .plot)) %>%
    dplyr::pull(plt)

ggsave('gradient_epistasis-s5.pdf', .plt, width=7, height=9.9)
