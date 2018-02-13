library(tidyverse)
library(wtl)
loadNamespace('cowplot')
#########1#########2#########3#########4#########5#########6#########7#########
source('~/git/likeligrid/R/common.R')

.datadir = '~/working/likeligrid/output'
.datadir = '~/working/likeligrid/output-pleiotropy'
.datadir = '~/working/likeligrid/output1115'
.datadir = '~/working/likeligrid/output1212'

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

read_top_noepi_gather = function(.x) {
  .result = list.files(.x, '\\.tsv\\.gz$', full.names=TRUE) %>% read_likeligrid()
  .epista_col = names(.result) %>% str_subset(':')
  if (length(.epista_col) == 1L) {
    .result = dplyr::rename(.result, epistasis = !!.epista_col)
  } else {
    .result = dplyr::mutate(epistasis = 1.0)
  }
  .top = .result %>% dplyr::top_n(1L, loglik) %>% dplyr::mutate(is_top=TRUE)
  .noepi = .result %>% dplyr::filter(epistasis == 1.0) %>% dplyr::top_n(1L, loglik) %>% dplyr::mutate(is_top=FALSE)
  .output = dplyr::bind_rows(.top, .noepi)
  if (!'pleiotropy' %in% names(.top)) {
    .output = .output %>% dplyr::mutate(pleiotropy = 1.0)
  }
  .output %>% tidyr::gather(pathway, value, -loglik, -is_top)
}
.metadata$indir[[2]] %>% read_top_noepi_gather()

.results_g5 = .metadata %>%
  dplyr::filter(s == 5L, gradient) %>%
  print() %>%
  dplyr::mutate(data= purrr::map(indir, read_top_noepi_gather)) %>%
  dplyr::select(-indir) %>%
  print()

.nested = .results_g5 %>%
  dplyr::select(-label) %>%
  tidyr::nest(-definition, -s, -gradient, -pleiotropy, -type) %>%
  dplyr::mutate(data = purrr::map(data, tidyr::unnest)) %>%
  print()

find_starting_point = function(x = .nested$data[[3]]) {x %>%
  dplyr::filter(!is_top) %>%
  tidyr::nest(-epistasis, -loglik, -is_top) %>%
  dplyr::top_n(1L, loglik) %>%
  dplyr::mutate(epistasis = NA_character_) %>%
  tidyr::unnest()
}

make_plottable = function(x = .nested$data[[3]]) {
  .start = x %>% find_starting_point()
  unique(.top$epistasis) %>% purrr::map_dfr(~{
    .start %>% dplyr::mutate(epistasis = .x)
  }) %>%
  dplyr::bind_rows(dplyr::filter(x, is_top))
}
make_plottable()

.plot = function(.tbl= make_plottable(.nested$data[[1]]), .title='') {
  message(.title)
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
