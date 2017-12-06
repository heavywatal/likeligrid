library(tidyverse)
library(wtl)
loadNamespace('cowplot')
#########1#########2#########3#########4#########5#########6#########7#########
source('~/git/likeligrid/R/common.R')

.datadir = '~/working/likeligrid/output'
.datadir = '~/working/likeligrid/output-pleiotropy'
.datadir = '~/working/likeligrid/output1115'

.metadata = read_metadata(.datadir) %>% print()

# #######1#########2#########3#########4#########5#########6#########7#########

.results = .metadata %>%
    dplyr::filter(is.na(epistasis)) %>%
    dplyr::filter(s == 4L) %>%
    # head(3) %>%
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
# summarize likelihood elevation in gradient mode

.results_g5 = .metadata %>%
  dplyr::filter(s == 5L, gradient) %>%
  print() %>%
  dplyr::mutate(data= purrr::map(indir, ~{
    .top = list.files(.x, '\\.tsv\\.gz$', full.names=TRUE) %>%
      read_likeligrid() %>%
      dplyr::top_n(1L, loglik)
    if (!any(str_detect(names(.top), ':'))) {
      .top = .top %>% dplyr::mutate(epistasis= 1.0)
    }
    if (!'pleiotropy' %in% names(.top)) {
      .top = .top %>% dplyr::mutate(pleiotropy= 1.0)
    }
    .top %>% tidyr::gather(pathway, value, -loglik)
  })) %>%
  dplyr::select(-indir) %>%
  print()

.nested = .results_g5 %>%
  tidyr::nest(-definition, -s, -gradient, -type) %>%
  print()

make_plottable = function(.typedata= .nested$data[[1]]) {
  .no_epi = {.typedata %>% dplyr::filter(is.na(epistasis))}$data[[1]]
  .typedata %>%
    dplyr::filter(!is.na(epistasis)) %>%
    dplyr::mutate(data= purrr::map_if(data, !pleiotropy, ~{dplyr::bind_rows(.no_epi, .x)})) %>%
    tidyr::unnest() %>%
    dplyr::mutate(pathway= str_replace(pathway, '\\w+:\\w+', 'epistasis'))
}
make_plottable()

.plot = function(.tbl= make_plottable(), .title='') {
  .start_loglik = min(.tbl$loglik)
  ggplot(.tbl, aes(value, loglik))+
  geom_hline(data=.const.chisq, aes(yintercept=.start_loglik + y, colour=alpha))+
  geom_path(size=1.3)+
  facet_grid(epistasis ~ pathway)+
  coord_cartesian(xlim=c(0, 2))+
  scale_x_continuous(breaks=c(0, 1, 2))+
  labs(title=.title)+
  theme_bw()+
  theme(legend.position='none')
}
.plot()

.plt = .nested %>%
  dplyr::mutate(data= purrr::map(data, make_plottable)) %>%
  dplyr::mutate(plt= purrr::map2(data, type, .plot)) %>%
  dplyr::pull(plt)

ggsave('gradient_epistasis-s5.pdf', .plt, width=7, height=9.9)
