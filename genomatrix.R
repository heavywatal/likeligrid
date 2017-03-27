library(pipeR)
library(stringr)
library(tidyverse)
library(jsonlite)
library(wtl)
loadNamespace('cowplot')
ggplot2::theme_set(wtl::theme_wtl())

repo = '~/git/innanlab'
cachedir = '~/Dropbox/working/innanlab'
raw_maf = file.path(cachedir, 'drivers.primary.PASS.extracted.tsv.gz') %>>% read_tsv() %>>%
    dplyr::select(-txid) %>>% (?.)

.ncases = raw_maf %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::count(cancer_type, sort=TRUE) %>>% (?.)

(.cancer_types = dplyr::filter(.ncases, nn > 400L)$cancer_type)

## TODO
# Assumption of impact?
# Count gene only once?
maf = raw_maf %>>%
    dplyr::filter(impact != 'LOW') %>>%
    dplyr::distinct(cancer_type, sample, symbol) %>>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>>%
    dplyr::arrange(cancer_type) %>>% (?.)

maf %>>% dplyr::filter(cancer_type %in% .cancer_types) %>>% dplyr::count(sample) %>>% (quantile(.$n, c(0.5, 0.95, 0.99, 1)))

(.hypermut_quantile = maf %>>% dplyr::count(sample) %>>% (quantile(.$n, c(0.5, 0.95, 0.99, 1))))
.hypermut_quantile['99%']
(.hypermut_df = tibble(p=names(.hypermut_quantile), x=.hypermut_quantile) %>>% dplyr::filter(p != '100%'))

.p = maf %>>%
    dplyr::mutate(cancer_type= factor(cancer_type, levels=.ncases$cancer_type)) %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::ungroup() %>>%
    ggplot(aes(n))+
    geom_histogram(bins=40)+
    geom_vline(data=.hypermut_df, aes(xintercept=x, colour=p))+
    scale_colour_manual(values=c('50%'='black', '95%'='brown', '99%'='tomato'))+
    facet_wrap(~cancer_type)+
    scale_x_log10(breaks=wtl::breaks_log10, labels=wtl::labels_log10)+
    labs(x='# mutated sites', y='# patients')+
    theme(legend.position=c(0.6, 0.05), legend.justification=c(0, 0), legend.direction='horizontal')
.p
ggsave('hist_mutated_sites.png', .p, width=12, height=7)

#########1#########2#########3#########4#########5#########6#########7#########

pathdefs = read_json('pathway-defs.json', simplifyVector=TRUE)
lengths(pathdefs)
pathdefs = pathdefs[c('vogelstein', 'vogelstein3', 'pereira', 'HALLMARK')]

.join = function(.genotypes, .pathdef) {
    .pd = map_df(.pathdef, ~tibble(symbol=.x), .id='pathway')
    dplyr::left_join(.genotypes, .pd, by='symbol')
}
.join(head(maf), pathdefs[['HALLMARK']])

join_pathdefs = function(.data) {
    purrr::map_df(pathdefs, .id='definition', ~.join(.data, .x)) %>>%
    dplyr::arrange(cancer_type, sample, symbol)
}
join_pathdefs(maf)

hypermutants = maf %>>%
  dplyr::count(cancer_type, sample) %>>%
  dplyr::filter(n > .hypermut_quantile['99%']) %>>%
  (unique(.$sample)) %>>% (?.)

major = maf %>>%
    dplyr::filter(cancer_type %in% .cancer_types) %>>%
    dplyr::filter(!sample %in% hypermutants) %>>%
    (?.)

shuffle = function(.data) {
    tidyr::nest(.data, -cancer_type) %>>%
    dplyr::mutate(data= purrr::map(data, ~{dplyr::mutate(.x, symbol= sample(symbol))})) %>>%
    tidyr::unnest()
}
major %>>% shuffle() %>>% dplyr::count(cancer_type, symbol)
major %>>% dplyr::count(cancer_type, symbol)

major %>>% dplyr::distinct(cancer_type, sample) %>>% dplyr::count(cancer_type)

#########1#########2#########3#########4#########5#########6#########7#########

.p = major %>>%
    # shuffle() %>>%
    join_pathdefs() %>>%
    # dplyr::filter(!is.na(pathway)) %>>%
    dplyr::count(definition, cancer_type, pathway) %>>%
    dplyr::ungroup() %>>%
    tidyr::nest(-definition) %>>%
    dplyr::mutate(plt= purrr::map2(data, definition, ~{
        ggplot(.x, aes(pathway, n))+
        geom_col()+
        facet_grid(~ cancer_type)+
        labs(title= .y)+
        wtl::theme_wtl(base_size=8)+
        theme(axis.ticks=element_blank(), panel.grid.major.x=element_blank())
    })) %>>%
    (plt) %>>% (cowplot::plot_grid(plotlist=., ncol=1))
.p
ggsave('pathway_freqs.png', .p, width=18, height=9)
# ggsave('pathway_freqs-narm.png', .p, width=18, height=9)
# ggsave('pathway_freqs-shuffled.png', .p, width=18, height=9)

.plot = function(.data, .title='', .max_g=max(.data$genes), .max_p=max(.data$pathways)) {
    .scatter = .data %>>%
        ggplot(aes(genes, pathways))+
        geom_jitter(width=0, height=0.4, alpha=0.25)+
        geom_vline(data=.hypermut_df, aes(xintercept=x, colour=p))+
        scale_colour_manual(values=c('50%'='black', '95%'='brown', '99%'='tomato'))+
        coord_cartesian(xlim=c(0, .max_g), ylim=c(0, .max_p))+
        theme(legend.position='none')
    .bar_g = .data %>>%
        ggplot(aes(genes))+
        geom_histogram(bins=30)+
        labs(title=.title)+
        coord_cartesian(xlim=c(0, .max_g))+
        theme(axis.title.x=element_blank())
    .bar_p = .data %>>%
        ggplot(aes(pathways))+
        geom_bar()+
        coord_flip(xlim=c(0, .max_p))+
        theme(axis.title.y=element_blank())
    cowplot::plot_grid(.bar_g, NULL, .scatter, .bar_p, rel_widths=c(2, 1), rel_heights=c(1, 2))
}

plot_mutdist = function(.data, .definition='HALLMARK', .quantile='99%', .head=12) {
    .max_g = .hypermut_quantile[.quantile]
    .data %>>%
    dplyr::filter(definition == .definition, genes <= .max_g) %>>%
    tidyr::nest(-cancer_type, -definition) %>>%
    head(.head) %>>%
    dplyr::mutate(plt= purrr::map2(data, paste(definition, cancer_type), .plot,
        .max_g=.max_g, .max_p=length(pathdefs[[.definition]]))) %>>%
    (cowplot::plot_grid(plotlist=.$plt, ncol=4))
}
# plot_mutdist(.counts.narm, .head=4)

.count_mutated = function(.data) {
    dplyr::full_join(by=c('cancer_type', 'sample'),
      .data %>>%
        dplyr::distinct(cancer_type, sample, symbol) %>>%
        dplyr::count(cancer_type, sample) %>>%
        dplyr::ungroup() %>>% dplyr::rename(genes=n),
      .data %>>%
        dplyr::count(cancer_type, sample, definition, pathway) %>>%
        dplyr::count(cancer_type, sample, definition) %>>%
        dplyr::ungroup() %>>% dplyr::rename(pathways=nn))
}

.counts = maf %>>% join_pathdefs() %>>% .count_mutated() %>>% (?.)
.counts.narm = maf %>>% join_pathdefs() %>>% dplyr::filter(!is.na(pathway)) %>>% .count_mutated(maf) %>>% (?.)

tidyr::crossing(def= names(pathdefs), p= c('95%', '99%', '100%')) %>>%
purrr::by_row(~{
    message(paste(.$def, .$p))
    plot_mutdist(.counts.narm, .$def, .$p)
}) %>>%
(ggsave('mutdist-narm.pdf', .$.out, width=16, height=10))


#########1#########2#########3#########4#########5#########6#########7#########
## Compare to Poisson and shuffled

.nest_by_definition = function(.data) {
    dplyr::filter(.data, !is.na(pathway)) %>>%
    dplyr::count(definition, cancer_type, sample, pathway) %>>%
    dplyr::ungroup() %>>%
    tidyr::nest(-definition) %>>% (?.)
}

.shuffled = major %>>%
    shuffle() %>>%
    join_pathdefs() %>>%
    .nest_by_definition() %>>% (?.)

.nested = major %>>%
    join_pathdefs() %>>%
    .nest_by_definition() %>>%
    dplyr::left_join(.shuffled %>>% dplyr::rename(shuffled=data),
      by=c('definition')) %>>%
    (?.)

## observed/shuffle vs poisson

make_poisson = function(.data) {
    nsam = dplyr::n_distinct(.data[['sample']]) #%>>% (?.)
    dplyr::count(.data, pathway, wt=n) %>>%
    dplyr::mutate(lambda= nn / nsam) %>>% #(?.) %>>%
    purrr::by_row(~{
        tibble(n= seq_len(40), y= dpois(n, .x$lambda) * nsam) %>>%
        dplyr::filter(y >= 1)
    }) %>>%
    dplyr::select(-nn, -lambda) %>>%
    unnest()
}
.nested$data[[1]] %>>% make_poisson() %>>% dplyr::summarise(muts= sum(n * y))
.nested$data[[1]] %>>% summarise(muts= sum(n))

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

.p = .nested %>>%
    purrr::by_row(~.plot(.x$data[[1]], .x$definition)) %>>%
    # purrr::by_row(~.plot(.x$shuffled[[1]], .x$definition)) %>>%
    (cowplot::plot_grid(plotlist=.$.out))
.p
ggsave('poisson_pathway-pancancer.pdf', .p, width=18, height=12)
# ggsave('poisson_pathway_shuffled-pancancer.pdf', .p, width=18, height=12)

.nested %>>% purrr::by_row(~{
    .outfile = sprintf('poisson_pathway-%s.pdf', .x$definition)
    # .outfile = sprintf('poisson_pathway_shuffled-%s.pdf', .x$definition)
    message(.outfile)
    .x$data[[1]] %>>%
    # .x$shuffled[[1]] %>>%
        tidyr::nest(-cancer_type) %>>%
        purrr::pmap(function(cancer_type, data) {
            .plot(data, cancer_type)
        }) %>>%
        (cowplot::plot_grid(plotlist=.)) %>>%
        (ggsave(.outfile, ., width=24, height=18))
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

.nested %>>%
    dplyr::rename(title= definition) %>>%
    purrr::pmap(.plot_obs_shuf) %>>%
    (cowplot::plot_grid(plotlist=.)) %>>%
    (ggsave('shuffle_pathway-pancancer.pdf', ., width=24, height=18))

.nested %>>% purrr::pmap(function(definition, data, shuffled) {
    .outfile = sprintf('shuffle_pathway-%s.pdf', definition)
    message(.outfile)
    .data = dplyr::left_join(by='cancer_type',
        data %>>% tidyr::nest(-cancer_type),
        shuffled %>>% tidyr::nest(-cancer_type) %>>% dplyr::rename(shuffled= data))
    .data %>>%
        dplyr::mutate(cancer_type= as.character(cancer_type)) %>>%
        dplyr::rename(title= cancer_type) %>>%
        purrr::pmap(.plot_obs_shuf)  %>>%
    (cowplot::plot_grid(plotlist=.)) %>>%
    (ggsave(.outfile, ., width=24, height=18))
})

#########1#########2#########3#########4#########5#########6#########7#########
# Summarize pathways and their mutations

.tidy = major %>>%
    join_pathdefs() %>>%
    dplyr::filter(!is.na(pathway)) %>>% (?.)

.mutcount = .tidy %>>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange(definition, cancer_type, pathway, desc(n)) %>>% (?.)

.summary_pathways = .mutcount %>>%
    dplyr::group_by(definition, cancer_type, pathway) %>>%
    dplyr::summarise(n_muts=sum(n), n_genes=n()) %>>%
    dplyr::arrange(definition, cancer_type, desc(n_muts)) %>>%
    dplyr::ungroup() %>>% (?.)

.usable_pathways = .summary_pathways %>>%
    dplyr::filter(n_muts >= 10, n_genes >=5) %>>%  #TODO: condition
    dplyr::group_by(definition, cancer_type) %>>%
    dplyr::top_n(8L, wt=n_muts) %>>%
    dplyr::ungroup() %>>% (?.)

.dirty = .mutcount %>>%
    dplyr::transmute(definition, cancer_type, pathway, value=sprintf('%s %4d', symbol, n)) %>>%
    dplyr::group_by(definition, cancer_type, pathway) %>>%
    dplyr::mutate(i=seq_along(pathway)) %>>%
    tidyr::spread(i, value) %>>%
    (dplyr::left_join(.summary_pathways, .)) %>>% (?.)

write_tsv(.dirty, 'summary-pathways.tsv', na='')

.dirty %>>%
    tidyr::nest(-definition) %>>%
    purrr::by_row(~{
        .outfile = sprintf('summary-pathways-%s.tsv', .$definition)
        message(.outfile)
        write_tsv(.$data[[1]], .outfile, na='')
    })

####

plot_pathwaypergene = function(.data) {
    .data %>>%
    dplyr::distinct(definition, cancer_type, pathway, symbol) %>>%
    dplyr::count(definition, cancer_type, symbol) %>>%
    dplyr::ungroup() %>>% (?.) %>>%
    ggplot(aes(n))+
    geom_bar()+
    facet_grid(definition ~ cancer_type)+
    labs(x='# pathways')+
    scale_x_continuous(breaks=c(5L, 10L))+
    theme(axis.ticks=element_blank())
}

.p = .tidy %>>% plot_pathwaypergene()
.p
ggsave('mutcount-pathways.pdf', .p, width=10, height=7)

drivers = file.path(repo, 'driver_genes.tsv') %>>% read_tsv() %>>% (?.)
known_drivers = drivers %>>% dplyr::filter(!is.na(role)) %>>% (?.)
cgc_drivers = drivers %>>% dplyr::filter(!is.na(cgc_vartype)) %>>% (?.)

.mutcount %>>%
    dplyr::filter(symbol %in% known_drivers$symbol) %>>%
    plot_pathwaypergene()

.mutcount %>>%
    dplyr::filter(symbol %in% cgc_drivers$symbol) %>>%
    plot_pathwaypergene()

.mutcount_major = dplyr::left_join(
    .usable_pathways %>>% tidyr::nest(-definition, -cancer_type, .key=usable),
    .mutcount %>>% tidyr::nest(-definition, -cancer_type, .key=mutcount)) %>>%
    dplyr::mutate(data= purrr::map2(usable, mutcount, ~{
        dplyr::filter(.y, pathway %in% .x$pathway)
    })) %>>%
    tidyr::unnest(data) %>>% (?.)

.p = .mutcount_major %>>% plot_pathwaypergene()
.p
ggsave('mutcount-pathways-major.pdf', .p, width=10, height=7)

.mutcount_major %>>%
dplyr::count(definition, cancer_type, symbol) %>>%
dplyr::filter(nn > 1L) %>>%
dplyr::arrange(definition, cancer_type, desc(nn)) %>>%
# dplyr::filter(definition == 'vogelstein') %>>%
dplyr::group_by(definition, cancer_type) %>>%
dplyr::mutate(i = seq_along(symbol)) %>>%
dplyr::ungroup() %>>%
dplyr::transmute(definition, cancer_type, i, value = sprintf('%s %2d', symbol, nn)) %>>%
tidyr::spread(i, value) %>>% (?.) %>>%
write_tsv('multi-pathways.tsv', na='')

#########1#########2#########3#########4#########5#########6#########7#########
## Generate genotype matrix

.matrices = .tidy %>>%
    dplyr::count(definition, cancer_type, sample, pathway) %>>%
    dplyr::ungroup() %>>%
    tidyr::nest(-definition, -cancer_type) %>>%
    dplyr::mutate(data= purrr::map(data, ~{
        .major_pathways = dplyr::count(.x, pathway, wt=n) %>>%
            dplyr::top_n(8L, wt=nn) %>>%
            (pathway)
        .x %>>%
        dplyr::filter(pathway %in% .major_pathways) %>>%
        tidyr::spread(pathway, n, fill=0L) %>>%
        dplyr::select(-sample)
    })) %>>% (?.)

.outdir = '~/Dropbox/working/likeligrid/pathtypes'
.matrices %>>% purrr::by_row(~{
    .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.tsv.gz', .$cancer_type, .$definition))
    message(.outfile)
    write_tsv(.$data[[1]], .outfile, na='')
})

# preserving all pathways
.matrices = major %>>%
    join_pathdefs() %>>%
    dplyr::filter(!is.na(pathway)) %>>%
    dplyr::count(definition, cancer_type, sample, pathway) %>>%
    dplyr::ungroup() %>>%
    tidyr::nest(-definition) %>>%
    dplyr::mutate(data= purrr::map(data, ~{
        .x %>>%
        tidyr::spread(pathway, n, fill=0L) %>>%
        dplyr::select(-sample) %>>%
        tidyr::nest(-cancer_type)
    })) %>>%
    tidyr::unnest() %>>% (?.)

.outdir = 'full-genotypes'
.matrices %>>% purrr::by_row(~{
    .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.tsv.gz', .$cancer_type, .$definition))
    message(.outfile)
    write_tsv(.$data[[1]], .outfile, na='')
})


#########1#########2#########3#########4#########5#########6#########7#########
## Make small data for exact model

.mutcount = .tidy %>>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>>%
    dplyr::filter(cumsum(n) / sum(n) < 0.9) %>>%
    dplyr::filter(n > 1) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange(definition, cancer_type, pathway, desc(n)) %>>% (?.)

.usable_pathways = .mutcount %>>%
    dplyr::group_by(definition, cancer_type, pathway) %>>%
    dplyr::summarise(n_muts=sum(n), n_genes=n()) %>>%
    dplyr::arrange(definition, cancer_type, desc(n_muts)) %>>% (?.) %>>%
    dplyr::filter(cumsum(n_muts) / sum(n_muts) < 0.9) %>>%
    # dplyr::filter(n_genes >=10) %>>%  #TODO: condition
    dplyr::group_by(definition, cancer_type) %>>%
    dplyr::top_n(4L, wt=n_muts) %>>%
    dplyr::ungroup() %>>% (?.)

.tiny = .usable_pathways %>>%
    dplyr::select(definition, cancer_type, pathway) %>>%
    dplyr::left_join(.mutcount %>>% dplyr::select(-n)) %>>%
    dplyr::left_join(.tidy) %>>%
    dplyr::arrange(definition, cancer_type, sample) %>>%
    (?.)

.tiny %>>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>>%
    dplyr::summarise(n_muts=sum(n), n_genes=n())

.p = .tiny %>>% plot_pathwaypergene()
.p
ggsave('mutcount-pathways-tiny.pdf', .p, width=10, height=7)

.dirty = .tiny %>>%
    dplyr::count(definition, cancer_type, pathway, symbol) %>>%
    dplyr::ungroup() %>>%
    dplyr::arrange(definition, cancer_type, pathway, desc(n)) %>>% (?.) %>>%
    dplyr::transmute(definition, cancer_type, pathway, value=sprintf('%s %4d', symbol, n)) %>>%
    dplyr::group_by(definition, cancer_type, pathway) %>>%
    dplyr::mutate(i=seq_along(pathway)) %>>%
    tidyr::spread(i, value) %>>%
    (dplyr::left_join(.usable_pathways, .)) %>>% (?.)
write_tsv(.dirty, 'usable-pathways-tiny.tsv', na='')

.outdir = '~/Dropbox/working/likeligrid/genotypes'
.tiny %>>%
    dplyr::distinct(definition, cancer_type, sample, symbol) %>>%
    dplyr::mutate(value = 1L) %>>%
    tidyr::nest(-definition, -cancer_type) %>>%
    # dplyr::filter(definition == 'vogelstein') %>>%
    # head(4) %>>%
    purrr::by_row(~{
        .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.tsv.gz', .$cancer_type, .$definition))
        tidyr::spread(.x$data[[1]], symbol, value, 0L) %>>%
        dplyr::select(-sample) %>>% #(?.)
        write_tsv(.outfile, na='')
    })

genotype2bits = function(.data) {
    .data %>>%
    dplyr::mutate(value=1L) %>>%
    tidyr::spread(symbol, value, 0L) %>>%
    tidyr::unite_('bits', names(.)[-1L], sep='')
}

.outdir = '~/Dropbox/working/likeligrid/genotypes'
.tiny %>>%
    tidyr::nest(-definition, -cancer_type) %>>%
    # head(2) %>>%
    purrr::by_row(~{
        .outfile = file.path(.outdir, sprintf('TCGA-%s-%s.json.gz', .$cancer_type, .$definition))
        message(.outfile)
        .data = .x$data[[1]] %>>% dplyr::mutate(symbol=as.factor(symbol))
        .path = .data %>>%
            dplyr::distinct(pathway, symbol) %>>%
            genotype2bits() #%>>% (?.)
        .sample = .data %>>%
            dplyr::distinct(sample, symbol) %>>%
            genotype2bits() #%>>% (?.)
        .con = gzfile(.outfile, 'w')
        list(symbol= levels(.data$symbol),
             pathway= .path$pathway,
             annotation= .path$bits,
             sample= .sample$bits) %>>%
        jsonlite::write_json(.con, pretty=TRUE)
        close(.con)
    })
