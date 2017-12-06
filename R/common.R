read_metadata = function(.dir) {
    tibble(indir= list.dirs(.dir, recursive=FALSE)) %>%
    dplyr::mutate(label= basename(indir)) %>%
    tidyr::separate(label, c('TCGA', 'type', 'definition', 'rest'), '-', remove=FALSE, extra='merge', fill='right') %>%
    dplyr::mutate(rest= paste0('-', rest),
      s= str_extract(rest, '(?<=-s)\\d') %>% parse_integer(),
      gradient= str_detect(rest, '-g'),
      epistasis= str_extract(rest, '(?<=-e)[^-]+'),
      pleiotropy= str_detect(rest, '-p')) %>%
    dplyr::select(-TCGA, -rest) %>%
    dplyr::filter(!is.na(type))
}

read_likeligrid = function(file) {
    read_tsv(file, comment='#', col_types=cols(.default='d'))
}

read_uniaxis = function(file) {
    read_likeligrid(file) %>%
    dplyr::select_if(function(.x) {dplyr::n_distinct(.x) > 1 || .x[1] < 0}) %>%
    tidyr::gather(pathway, value, -loglik)
}

read_limit_max = function(file) {
    .pathway = str_extract(file, '(?<=limit-)[^.]+(?=_[LU]\\.tsv)')
    read_likeligrid(file) %>%
    dplyr::select(loglik, one_of(.pathway)) %>%
    setNames(c('loglik', 'value')) %>%
    dplyr::mutate(pathway=.pathway) %>%
    dplyr::group_by(value) %>%
    dplyr::top_n(1, wt=loglik) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(value)
}

read_all_uniaxis = function(directory) {
    list.files(directory, 'uniaxis.+\\.gz$', full.names=TRUE) %>%
    purrr::map_dfr(read_uniaxis)
}

read_all_limit = function(directory) {
    list.files(directory, 'limit.+\\.gz$', full.names=TRUE) %>%
    purrr::map_dfr(read_limit_max)
}

# Put interaction term at last of facet_wrap
tilde_epistasis = function(.names) {
    str_replace(.names, '(\\w+:\\w+|pleiotropy)', '~\\1')
}
# c('A', 'B', 'A:B', 'pleiotropy') %>% tilde_epistasis()

.const.chisq = tibble(alpha=c(0.95, 0.99), y=qchisq(alpha, 1) * 0.5)

plot_landscape = function(uniaxis, limits=tibble(), gradient=tibble(), label='', ...) {
    if (gradient %>% inherits('list')) {
        return(purrr::map(gradient, ~plot_landscape(uniaxis, limits, .x, label, ...)))
    }
    .mle = dplyr::top_n(uniaxis, 1L, loglik) %>% dplyr::distinct(loglik, pathway, .keep_all=TRUE)
    .max_loglik = .mle$loglik[1]
    .max_y = .max_loglik
    .p = uniaxis %>%
        mutate(pathway= tilde_epistasis(pathway)) %>%
        ggplot(aes(value, loglik))+
        geom_line(colour='#66cc66')
    if (nrow(limits) > 0L) {
        .p = .p + geom_point(data=limits %>% mutate(pathway= tilde_epistasis(pathway)), alpha=0.5)
    }
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
          mutate(pathway= tilde_epistasis(pathway)) %>%
          dplyr::mutate(loglik=loglik+.offset)
        .p = .p + geom_point(data=.data_g, colour='#33cc99', alpha=0.5)
        .max_y = max(.data_g$loglik)
    }
    .p+
    geom_vline(data=.mle, aes(xintercept=value), colour='salmon')+
    geom_hline(data=.const.chisq, aes(yintercept=.max_loglik - y, colour=alpha))+
    coord_cartesian(ylim=c(.max_y, .max_loglik - qchisq(0.95, 1)))+
    scale_x_continuous(breaks=c(0, 1, 2), limits=c(0, 2), expand=c(0, 0))+
    labs(title=label)+
    facet_wrap(~pathway)+
    wtl::theme_wtl()+
    theme(legend.position='none', panel.grid.major.y=element_blank(),
      panel.grid.major.x=element_line(colour='#aaaaaa'),
      axis.ticks=element_blank(), axis.title.x=element_blank())
}
