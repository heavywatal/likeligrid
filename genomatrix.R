library(stringr)
library(tidyverse)
library(jsonlite)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())

repo = '~/git/innanlab'
cachedir = '~/Dropbox/working/innanlab'
raw_maf = file.path(cachedir, 'drivers.primary.PASS.extracted.tsv.gz') %>% read_tsv() %>%
    dplyr::select(-txid) %>% print()

.ncases = raw_maf %>%
    dplyr::count(cancer_type, sample) %>%
    dplyr::count(cancer_type, sort=TRUE) %>% print()

(.cancer_types400 = dplyr::filter(.ncases, nn > 400L)$cancer_type)

# Kandoth 2013 Nature, Tamborero 2013 Sci Rep
major_cancer_types = c('BRCA', 'LUAD', 'LUSC', 'UCEC', 'GBM', 'HNSC', 'COAD', 'READ', 'BLCA', 'KIRC', 'OV', 'LAML')

.ncases %>% dplyr::filter(cancer_type %in% major_cancer_types)
union(.cancer_types400, major_cancer_types)
intersect(.cancer_types400, major_cancer_types)

(.cancer_types = str_replace(major_cancer_types, 'LAML', 'STAD'))

## TODO
# Assumption of impact?
# Count gene only once?
maf = raw_maf %>%
    dplyr::filter(cancer_type %in% .cancer_types) %>%
    dplyr::filter(impact != 'LOW') %>%
    dplyr::distinct(cancer_type, sample, symbol) %>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>%
    dplyr::arrange(cancer_type) %>% print()

(.hypermut_quantile = maf %>% dplyr::count(sample) %>% {quantile(.$n, c(0.5, 0.95, 0.99, 1))})
.hypermut_quantile['99%']
(.hypermut_df = tibble(p=names(.hypermut_quantile), x=.hypermut_quantile) %>% dplyr::filter(p != '100%'))

#########1#########2#########3#########4
## Filter by annotation?

drivers = file.path(repo, 'driver_genes.tsv') %>% read_tsv() %>% print()
known_drivers = drivers %>% dplyr::filter(!is.na(role)) %>% print()
cgc_drivers = drivers %>% dplyr::filter(!is.na(cgc_vartype)) %>% print()

#########1#########2#########3#########4#########5#########6#########7#########

pathdefs = read_json('pathway-defs.json', simplifyVector=TRUE)
lengths(pathdefs)
pathdefs = pathdefs[c('vogelstein', 'vogelstein3', 'pereira', 'HALLMARK')]

pathdefs_df = pathdefs %>%
    purrr::map(tibble::enframe, 'pathway', 'symbol') %>%
    tibble::enframe('definition') %>%
    tidyr::unnest() %>%
    dplyr::mutate(symbol= purrr::map(symbol, ~{.x[.x %in% cgc_drivers$symbol]})) %>%
    dplyr::filter(purrr::map_lgl(symbol, ~{length(.x) > 1})) %>%
    print()

pathdefs_cancer = pathdefs_df %>%
    tidyr::nest(-definition) %>%
    dplyr::mutate(data= purrr::map(data, tibble::deframe)) %>%
    tibble::deframe() %>%
    print()

pathdefs %>% purrr::map(~{flatten_chr(.x) %>% unique()}) %>% lengths()
pathdefs_cancer %>% purrr::map(~{flatten_chr(.x) %>% unique()}) %>% lengths()

.join = function(.genotypes, .pathdef) {
    .pd = map_df(.pathdef, ~tibble(symbol=.x), .id='pathway')
    dplyr::left_join(.genotypes, .pd, by='symbol')
}
.join(head(maf), pathdefs[['HALLMARK']])

join_pathdefs = function(.data) {
    purrr::map_df(pathdefs_cancer, .id='definition', ~.join(.data, .x)) %>%
    dplyr::arrange(cancer_type, sample, symbol)
}
join_pathdefs(maf)

#########1#########2#########3#########4#########5#########6#########7#########
# Summarize pathways and their mutations

show_ngene_per_def = function(.data) {.data %>%
    dplyr::distinct(cancer_type, definition, symbol) %>%
    dplyr::count(cancer_type, definition) %>%
    dplyr::arrange(definition, cancer_type) %>%
    wtl::max_print()
}

hypermutants = maf %>%
  dplyr::count(cancer_type, sample) %>%
  dplyr::filter(n > .hypermut_quantile['99%']) %>%
  {unique(.$sample)} %>% print()

hypermutants_each = maf %>%
    dplyr::count(cancer_type, sample) %>%
    dplyr::arrange(cancer_type, n) %>%
    dplyr::mutate(i = seq_len(n())) %>%
    dplyr::filter(i / n() > 0.90) %>%
    {.$sample} %>% print()

.tidy = maf %>%
    dplyr::filter(!sample %in% hypermutants_each) %>%
    ################################################## exclude hypermutants
    join_pathdefs() %>%
    dplyr::filter(!is.na(pathway)) %>% print()

.tidy %>% show_ngene_per_def()

.mutcount = .tidy %>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>%
    dplyr::arrange(definition, cancer_type, pathway, desc(n)) %>%
    ############################################# exclude less-mutated genes
    dplyr::filter(cumsum(n) / sum(n) < 0.95) %>%
    dplyr::filter(n > 1) %>%
    dplyr::ungroup() %>% print()

.mutcount %>% show_ngene_per_def()

.summary_pathways = .mutcount %>%
    dplyr::group_by(definition, cancer_type, pathway) %>%
    dplyr::summarise(n_muts=sum(n), n_genes=n()) %>%
    dplyr::arrange(definition, cancer_type, desc(n_muts)) %>%
    dplyr::ungroup() %>% print()

.summary_pathways %>%
    dplyr::filter(definition == 'vogelstein') %>%
    dplyr::group_by(definition, cancer_type) %>%
    dplyr::mutate(i= seq_len(n()), csum= cumsum(n_muts) / sum(n_muts)) %>%
    wtl::less()

.major_pathways = .summary_pathways %>%
    dplyr::group_by(definition, cancer_type) %>%
    dplyr::arrange(definition, cancer_type, desc(n_muts)) %>%
    dplyr::mutate(i= seq_len(n())) %>%
    ################################### exclude less-mutated pathways
    dplyr::filter(cumsum(n_muts) / sum(n_muts) < 0.95 | i < 4L) %>%
    dplyr::filter(n_genes > 1L) %>%  #TODO: condition
    dplyr::select(-i) %>%
    dplyr::ungroup() %>% print()

wtl::less(.major_pathways)

.major_pathways %>% dplyr::filter(n_genes < 3)

.major_pathways_tiny = .major_pathways %>%
    dplyr::group_by(definition, cancer_type) %>%
    dplyr::top_n(4L, wt=n_muts) %>%
    dplyr::ungroup() %>% print()

.major_pathways_manual = .major_pathways %>%
    dplyr::group_by(definition, cancer_type) %>%
    dplyr::top_n(6L, wt=n_muts) %>%
    dplyr::ungroup() %>% print()

.major_pathways_manual %>%
    dplyr::group_by(definition, cancer_type) %>%
    dplyr::filter(n() > 6)

.major = .major_pathways %>%
    dplyr::select(definition, cancer_type, pathway) %>%
    dplyr::left_join(.mutcount %>% dplyr::select(-n)) %>%
    dplyr::left_join(.tidy) %>%
    dplyr::arrange(definition, cancer_type, sample) %>%
    print()

.tiny = .major_pathways_tiny %>%
    dplyr::select(definition, cancer_type, pathway) %>%
    dplyr::left_join(.mutcount %>% dplyr::select(-n)) %>%
    dplyr::left_join(.tidy) %>%
    dplyr::arrange(definition, cancer_type, sample) %>%
    print()

.manual = .major_pathways_manual %>%
    dplyr::select(definition, cancer_type, pathway) %>%
    dplyr::left_join(.mutcount %>% dplyr::select(-n)) %>%
    dplyr::left_join(.tidy) %>%
    dplyr::arrange(definition, cancer_type, sample) %>%
    print()

.major %>% show_ngene_per_def() %>% write_tsv('ngene_per_def_major.tsv')
.tiny %>% show_ngene_per_def() %>% write_tsv('ngene_per_def_tiny.tsv')
.manual %>% show_ngene_per_def() %>% write_tsv('ngene_per_def_manual.tsv')

summarize_genes = function(.data, .def=c('vogelstein')) {.data %>%
    dplyr::filter(definition %in% .def) %>%
    dplyr::group_by(cancer_type, definition, symbol) %>%
    dplyr::summarize(
        nmut= n_distinct(sample),
        npath= n_distinct(pathway),
        pathways= str_c(unique(pathway), collapse=':')) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(cancer_type, definition, desc(nmut), symbol)
}
summarize_genes(.major) %>% less()
summarize_genes(.tiny) %>% less()
summarize_genes(.major, 'HALLMARK')
summarize_genes(.manual) %>% less()

.summary = summarize_genes(.manual) %>%
    dplyr::select(-definition) %>%
    dplyr::group_by(cancer_type) %>%
    dplyr::top_n(20L, nmut) %>%
    print() %>%
    write_tsv('summary-genes.tsv', na='')

loadNamespace('gridExtra')
.summary %>%
    tidyr::nest(-cancer_type) %>%# head(4) %>%
    dplyr::mutate(grob= purrr::map2(cancer_type, data, ~{
        .title = grid::textGrob(.x)
        .table = gridExtra::tableGrob(.y, theme=gridExtra::ttheme_minimal())
        .table = gtable::gtable_add_rows(.table, grid::unit(2, 'lines'), 0)
        gtable::gtable_add_grob(.table, .title, 1, 2, clip='off')
    })) %>%
    {gridExtra::marrangeGrob(.$grob, top=1, nrow=1, ncol=1)} %>%
    {ggsave('summary-genes.pdf', ., width=4, height=7)}

.dirty = .mutcount %>%
    dplyr::transmute(definition, cancer_type, pathway, value=sprintf('%s %4d', symbol, n)) %>%
    dplyr::group_by(definition, cancer_type, pathway) %>%
    dplyr::mutate(i=seq_along(pathway)) %>%
    tidyr::spread(i, value) %>%
    {dplyr::left_join(.summary_pathways, .)} %>% print()
write_tsv(.dirty, 'summary-pathways.tsv', na='')

.dirty %>%
    tidyr::nest(-definition) %>%
    purrr::pwalk(function(definition, data) {
        .outfile = sprintf('summary-pathways-%s.tsv', definition)
        message(.outfile)
        write_tsv(data, .outfile, na='')
    })

.dirty = .manual %>%
# .dirty = .major %>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(definition, cancer_type, pathway, desc(n)) %>% print() %>%
    dplyr::transmute(definition, cancer_type, pathway, value=sprintf('%s %4d', symbol, n)) %>%
    dplyr::group_by(definition, cancer_type, pathway) %>%
    dplyr::mutate(i=seq_along(pathway)) %>%
    tidyr::spread(i, value) %>%
    {dplyr::left_join(.major_pathways, .)} %>% print()
# write_tsv(.dirty, 'summary-pathways-major.tsv', na='')
write_tsv(.dirty, 'summary-pathways-manual.tsv', na='')

####

hist_pleiotropy = function(.data) {
    .data %>%
    dplyr::distinct(definition, cancer_type, pathway, symbol) %>%
    dplyr::count(definition, cancer_type, symbol) %>%
    dplyr::ungroup() %>% print() %>%
    ggplot(aes(n))+
    geom_bar()+
    facet_grid(definition ~ cancer_type)+
    labs(x='# pathways')+
    scale_x_continuous(breaks=seq.int(2L, 10L, 2L))+
    theme(axis.ticks=element_blank())
}

.p = .major %>% hist_pleiotropy()
ggsave('hist_pleiotropy-major.pdf', .p, width=10, height=7)

.p = .tiny %>% hist_pleiotropy()
ggsave('hist_pleiotropy-tiny.pdf', .p, width=10, height=7)

.p = .manual %>% hist_pleiotropy()
ggsave('hist_pleiotropy-manual.pdf', .p, width=10, height=7)

.mutcount %>%
dplyr::count(definition, cancer_type, symbol) %>%
dplyr::filter(nn > 1L) %>%
dplyr::arrange(definition, cancer_type, desc(nn)) %>%
# dplyr::filter(definition == 'vogelstein') %>%
dplyr::group_by(definition, cancer_type) %>%
dplyr::mutate(i = seq_along(symbol)) %>%
dplyr::ungroup() %>%
dplyr::transmute(definition, cancer_type, i, value = sprintf('%s %2d', symbol, nn)) %>%
tidyr::spread(i, value) %>% print() %>%
write_tsv('pleiotropic-genes.tsv', na='')

#########1#########2#########3#########4#########5#########6#########7#########
## Make binary genotypes

genotype2bits = function(.data) {
    .data %>%
    dplyr::mutate(value=1L) %>%
    tidyr::spread(symbol, value, 0L) %>%
    tidyr::unite_('bits', names(.)[-1L], sep='')
}

# .outdir = sprintf('~/Dropbox/working/likeligrid/genotypes/major-%s', wtl::today())
# .outdir = sprintf('~/Dropbox/working/likeligrid/genotypes/tiny-%s', wtl::today())
.outdir = sprintf('~/Dropbox/working/likeligrid/genotypes/manual-%s', wtl::today())
dir.create(.outdir, mode='0755')
.major %>%
# .tiny %>%
# .manual %>%
    tidyr::nest(-definition, -cancer_type) %>%
    dplyr::mutate(cancer_type= as.character(cancer_type)) %>%
    # head(2) %>%
    purrr::pwalk(function(definition, cancer_type, data) {
        .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.json.gz', cancer_type, definition))
        message(.outfile)
        .data = data %>% dplyr::mutate(symbol=as.factor(symbol))
        .path = .data %>%
            dplyr::distinct(pathway, symbol) %>%
            genotype2bits()
        .sample = .data %>%
            dplyr::distinct(sample, symbol) %>%
            genotype2bits()
        .con = gzfile(.outfile, 'w')
        list(symbol= levels(.data$symbol),
             pathway= .path$pathway,
             annotation= .path$bits,
             sample= .sample$bits) %>%
        jsonlite::write_json(.con, pretty=TRUE)
        close(.con)
    })

# #######1#########2#########3#########4#########5#########6#########7#########
## Visualize genotype matrix

.fun = function(definition, cancer_type, data) {
    .title = sprintf('TCGA-%s-%s', cancer_type, definition)
    message(.title)
    .data = data %>%
        dplyr::group_by(sample, symbol) %>%
        dplyr::summarize(pathway= paste0(pathway, collapse=':')) %>%
        dplyr::ungroup()
    .symbol = .data %>%
        dplyr::count(symbol) %>%
        dplyr::arrange(desc(n)) %>%
        {.$symbol}
    .sample = .data %>%
        dplyr::mutate(symbol= factor(symbol, levels=.symbol)) %>%
        dplyr::arrange(symbol) %>%
        dplyr::distinct(sample) %>%
        {.$sample}
    .data %>%
        dplyr::mutate(sample= factor(sample, levels=.sample)) %>%
        .plot(.title)
}

.plot = function(.data, .title) {
    ggplot(.data, aes(sample, symbol))+
    geom_tile()+
    # blur with geom_raster()?
    facet_grid(pathway ~ ., scales='free_y', space='free_y', switch='y')+
    scale_x_discrete(position='top')+
    labs(x=NULL, title=.title)+
    wtl::theme_wtl()+
    theme(panel.grid=element_blank(), axis.ticks=element_blank(),
        axis.text.x=element_blank(),
        strip.text.y=element_text(angle=180))
}

.manual %>%
    dplyr::distinct() %>%
    tidyr::nest(-definition, -cancer_type) %>%
    dplyr::mutate(cancer_type= as.character(cancer_type)) %>%
    dplyr::filter(definition == 'vogelstein') %>% #head(2) %>%
    purrr::pmap(.fun) %>%
    {ggsave('TCGA-genotype-matrix.pdf', ., width=18, height=10)}

#########1#########2#########3#########4#########5#########6#########7#########
## Visualize mutation distribution

.plot = function(.data, .title='', .max_g=max(.data$genes), .max_p=max(.data$pathways) + 0.5, .vline=FALSE) {
    .scatter = .data %>%
        ggplot(aes(genes, pathways))+
        geom_jitter(width=0, height=0.25, alpha=0.25)+
        coord_cartesian(xlim=c(0, .max_g), ylim=c(0, .max_p))+
        theme(legend.position='none')
    if (.vline) {
        .scatter = .scatter +
        geom_vline(data=.hypermut_df, aes(xintercept=x, colour=p))+
        scale_colour_manual(values=c('50%'='black', '95%'='brown', '99%'='tomato'))
    }
    .bar_g = .data %>%
        ggplot(aes(genes))+
        # geom_histogram(bins=30)+
        geom_bar()+
        labs(title=.title)+
        coord_cartesian(xlim=c(0, .max_g))+
        theme(axis.title.x=element_blank())
    .bar_p = .data %>%
        ggplot(aes(pathways))+
        geom_bar()+
        coord_flip(xlim=c(0, .max_p))+
        theme(axis.title.y=element_blank())
    cowplot::plot_grid(.bar_g, NULL, .scatter, .bar_p, rel_widths=c(2, 1), rel_heights=c(1, 2))
}

plot_mutdist = function(.data, .definition='HALLMARK', .quantile='99%', .head=12) {
    .max_g = .hypermut_quantile[.quantile]
    .data %>%
    dplyr::filter(definition == .definition, genes <= .max_g) %>%
    tidyr::nest(-cancer_type, -definition) %>%
    head(.head) %>%
    dplyr::mutate(plt= purrr::map2(data, paste(definition, cancer_type), .plot,
        .max_g=.max_g, .max_p=length(pathdefs[[.definition]]), .vline=TRUE)) %>%
    {cowplot::plot_grid(plotlist=.$plt, ncol=4)}
}
# plot_mutdist(.counts.narm, .head=4)

.count_mutated = function(.data) {
    dplyr::full_join(by=c('cancer_type', 'sample'),
      .data %>%
        dplyr::distinct(cancer_type, sample, symbol) %>%
        dplyr::count(cancer_type, sample) %>%
        dplyr::ungroup() %>% dplyr::rename(genes=n),
      .data %>%
        dplyr::count(cancer_type, sample, definition, pathway) %>%
        dplyr::count(cancer_type, sample, definition) %>%
        dplyr::ungroup() %>% dplyr::rename(pathways=nn))
}

.counts = maf %>% join_pathdefs() %>% .count_mutated() %>% print()
.counts.narm = maf %>% join_pathdefs() %>% dplyr::filter(!is.na(pathway)) %>% .count_mutated() %>% print()

tidyr::crossing(def= names(pathdefs), p= c('95%', '99%', '100%')) %>% head(2) %>%
purrr::pmap(function(def, p) {
    message(paste(def, p))
    plot_mutdist(.counts.narm, def, p)
}) %>%
{ggsave('mutdist-narm.pdf', ., width=16, height=10)}

# .major %>%
# .tiny %>%
.manual %>%
.count_mutated() %>%
tidyr::nest(-cancer_type, -definition) %>%
dplyr::mutate(plt= purrr::map2(data, paste(definition, cancer_type), .plot)) %>%
dplyr::arrange(definition, cancer_type) %>%
dplyr::group_by(definition) %>%
purrr::by_slice(~{
    cowplot::plot_grid(plotlist=.$plt, ncol=4)
}) %>%
{ggsave('mutdist-manual.pdf', .$.out, width=16, height=10)}
# {ggsave('mutdist-tiny.pdf', .$.out, width=16, height=10)}
# {ggsave('mutdist-major.pdf', .$.out, width=16, height=10)}

.major %>%
.manual %>%
    dplyr::count(cancer_type, definition, sample, pathway) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-cancer_type, -definition) %>%
    dplyr::arrange(definition, cancer_type) %>%
    purrr::pmap(function(definition, cancer_type, data) {
        .title = paste(cancer_type, definition)
        message(.title)
        ggplot(data, aes(n))+
        geom_bar()+
        facet_wrap(~pathway)+
        labs(title=.title)+
        wtl::theme_wtl()
    }) %>%
    {ggsave('hist-s-pathway-manual.pdf', ., width=12, height=12)}
    # {ggsave('hist-s-pathway-major.pdf', ., width=12, height=12)}

#########1#########2#########3#########4#########5#########6#########7#########
## Shuffle

nonhypermut = maf %>%
    dplyr::filter(cancer_type %in% .cancer_types) %>%
    dplyr::filter(!sample %in% hypermutants) %>%
    print()

shuffle = function(.data) {
    tidyr::nest(.data, -cancer_type) %>%
    dplyr::mutate(data= purrr::map(data, ~{dplyr::mutate(.x, symbol= sample(symbol))})) %>%
    tidyr::unnest()
}
nonhypermut %>% shuffle() %>% dplyr::count(cancer_type, symbol)
nonhypermut %>% dplyr::count(cancer_type, symbol)

nonhypermut %>% dplyr::distinct(cancer_type, sample) %>% dplyr::count(cancer_type)

.p = nonhypermut %>%
    # shuffle() %>%
    join_pathdefs() %>%
    # dplyr::filter(!is.na(pathway)) %>%
    dplyr::count(definition, cancer_type, pathway) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-definition) %>%
    dplyr::mutate(plt= purrr::map2(data, definition, ~{
        ggplot(.x, aes(pathway, n))+
        geom_col()+
        facet_grid(~ cancer_type)+
        labs(title= .y)+
        wtl::theme_wtl(base_size=8)+
        theme(axis.ticks=element_blank(), panel.grid.major.x=element_blank())
    })) %>%
    {.$plt} %>% {cowplot::plot_grid(plotlist=., ncol=1)}
.p
ggsave('pathway_freqs.png', .p, width=18, height=9)
# ggsave('pathway_freqs-narm.png', .p, width=18, height=9)
# ggsave('pathway_freqs-shuffled.png', .p, width=18, height=9)

#########1#########2#########3#########4#########5#########6#########7#########
## Compare to Poisson and shuffled

.nest_by_definition = function(.data) {
    dplyr::filter(.data, !is.na(pathway)) %>%
    dplyr::count(definition, cancer_type, sample, pathway) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-definition) %>% print()
}

.shuffled = nonhypermut %>%
    shuffle() %>%
    join_pathdefs() %>%
    .nest_by_definition() %>% print()

.nested = nonhypermut %>%
    join_pathdefs() %>%
    .nest_by_definition() %>%
    dplyr::left_join(.shuffled %>% dplyr::rename(shuffled=data),
      by=c('definition')) %>%
    print()

## observed/shuffle vs poisson

make_poisson = function(.data) {
    nsam = dplyr::n_distinct(.data[['sample']])
    dplyr::count(.data, pathway, wt=n) %>%
    dplyr::mutate(lambda= nn / nsam) %>%
    purrr::pmap_df(function(lambda, ...) {
        tibble(n= seq_len(40), y= dpois(n, lambda) * nsam) %>%
        dplyr::filter(y >= 1)
    }) %>%
    dplyr::select(-nn, -lambda) %>%
    unnest()
}
.nested$data[[1]] %>% make_poisson() %>% dplyr::summarise(muts= sum(n * y))
.nested$data[[1]] %>% summarise(muts= sum(n))

.plot = function(.data, .title) {
    .pois_df = make_poisson(.data)
    ggplot(.data, aes(n))+
    geom_bar(fill='dodgerblue', alpha=0.5)+
    geom_col(data=.pois_df, aes(n, y), fill='tomato', alpha=0.5)+
    facet_wrap(~ pathway)+
    labs(title= .title)+
    # scale_y_log10(breaks=wtl::breaks_log10, labels=wtl::labels_log10)+
    wtl::theme_wtl(base_size=10)+
    theme(axis.ticks=element_blank(), panel.grid.major.x=element_blank())
}
.plot(.nested$data[[3]], .nested$definition[[3]])

.p = .nested %>%
    purrr::pmap(function(definition, data, ...) {.plot(data, definition)}) %>%
    # purrr::pmap(function(definition, shuffled, ...) {.plot(shuffled, definition)}) %>%
    {cowplot::plot_grid(plotlist=.$.out)}
.p
ggsave('poisson_pathway-pancancer.pdf', .p, width=18, height=12)
# ggsave('poisson_pathway_shuffled-pancancer.pdf', .p, width=18, height=12)

.nested %>% purrr::pwalk(function(definition, data) {
    .outfile = sprintf('poisson_pathway-%s.pdf', definition)
    # .outfile = sprintf('poisson_pathway_shuffled-%s.pdf', .x$definition)
    message(.outfile)
    data %>%
    # .x$shuffled[[1]] %>%
        tidyr::nest(-cancer_type) %>%
        purrr::pmap(function(cancer_type, data) {
            .plot(data, cancer_type)
        }) %>%
        {cowplot::plot_grid(plotlist=.)} %>%
        {ggsave(.outfile, ., width=24, height=18)}
})

## observed vs shuffle

.plot_obs_shuf = function(title, data, shuffled) {
    message(title)
    ggplot(data, aes(n))+
    # geom_col(data=.pois_df, aes(n, y), fill='forestgreen', alpha=0.3)+
    geom_bar(fill='dodgerblue', alpha=0.4)+
    geom_bar(data=shuffled, fill='tomato', alpha=0.4)+
    facet_wrap(~ pathway)+
    labs(title= title)+
    # scale_y_log10(breaks=wtl::breaks_log10, labels=wtl::labels_log10)+
    wtl::theme_wtl(base_size=10)+
    theme(axis.ticks=element_blank(), panel.grid.major.x=element_blank())
}

.nested %>%
    dplyr::rename(title= definition) %>%
    purrr::pmap(.plot_obs_shuf) %>%
    {cowplot::plot_grid(plotlist=.)} %>%
    {ggsave('shuffle_pathway-pancancer.pdf', ., width=24, height=18)}

.nested %>% purrr::pmap(function(definition, data, shuffled) {
    .outfile = sprintf('shuffle_pathway-%s.pdf', definition)
    message(.outfile)
    .data = dplyr::left_join(by='cancer_type',
        data %>% tidyr::nest(-cancer_type),
        shuffled %>% tidyr::nest(-cancer_type) %>% dplyr::rename(shuffled= data))
    .data %>%
        dplyr::mutate(cancer_type= as.character(cancer_type)) %>%
        dplyr::rename(title= cancer_type) %>%
        purrr::pmap(.plot_obs_shuf)  %>%
    {cowplot::plot_grid(plotlist=.)} %>%
    {ggsave(.outfile, ., width=24, height=18)}
})

#########1#########2#########3#########4#########5#########6#########7#########
## Generate pathtype matrix

.matrices = .manual %>%
    dplyr::count(definition, cancer_type, sample, pathway) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-definition, -cancer_type) %>%
    dplyr::mutate(data= purrr::map(data, ~{
        .x %>%
        tidyr::spread(pathway, n, fill=0L) %>%
        dplyr::select(-sample)
    })) %>% print()

.outdir = '~/Dropbox/working/likeligrid/pathtypes'
.matrices %>% purrr::pwalk(function(definition, cancer_type, data) {
    .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.tsv.gz', cancer_type, definition))
    message(.outfile)
    write_tsv(data, .outfile, na='')
})

# preserving all pathways
.matrices = .tidy %>%
    dplyr::count(definition, cancer_type, sample, pathway) %>%
    dplyr::ungroup() %>%
    tidyr::nest(-definition) %>%
    dplyr::mutate(data= purrr::map(data, ~{
        .x %>%
        tidyr::spread(pathway, n, fill=0L) %>%
        dplyr::select(-sample) %>%
        tidyr::nest(-cancer_type)
    })) %>%
    tidyr::unnest() %>% print()

.outdir = 'full-pathtypes'
.matrices %>% purrr::pwalk(function(definition, cancer_type, data) {
    .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.tsv.gz', cancer_type, definition))
    message(.outfile)
    write_tsv(data, .outfile, na='')
})
