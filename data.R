#########1#########2#########3#########4#########5#########6#########7#########
## Pereira 2016 Nat Comm

pereiradir = '~/git/innanlab/pereira'
drivers = read_tsv(file.path(pereiradir, 'drivers40.tsv')) %>>% (?.)
datadir = file.path(pereiradir, 'mutationalProfiles/Data')
list.files(datadir)

patients = file.path(datadir, 'patientData.txt') %>>% read_tsv() %>>% (?.)

mutations_all = file.path(datadir, 'somaticMutations_incNC.txt') %>>% read_tsv() %>>% (?.)
mutations_all %>>% dplyr::count(location, mutationType)

mutations = file.path(datadir, 'somaticMutations.txt') %>>% read_tsv() %>>% (?.)
mutations %>>% dplyr::count(location, mutationType)

observed = mutations %>>%
    dplyr::left_join(drivers, by='gene') %>>%
    dplyr::transmute(sample, gene,
        pathway=coalesce(pathway, 'Other'),
        mutype=ifelse(str_detect(location ,'splicing'), 'truncating',
               ifelse(str_detect(mutationType, '^missense|^inframe'), 'missense', 'truncating'))) %>>%
    dplyr::left_join(patients %>>% dplyr::select(sample, erStatus, her2Status), by='sample') %>>% (?.)

observed %>>% dplyr::count(erStatus, her2Status)

genotype_pereira = observed %>>%
    dplyr::filter(erStatus == 'POS') %>>%
    dplyr::count(sample, pathway) %>>%
    tidyr::spread(pathway, n, 0L) %>>%
    dplyr::ungroup() %>>% dplyr::select(-sample) %>>%
    {ifelse(. > 0L, 1L, 0L)} %>>% as_tibble() %>>%
#    dplyr::mutate(epsilon=1L) %>>%
    (?.)

genotype_pereira %>>% summarise_all(mean)

.outdir = '~/Dropbox/working/cancer/genotypes'
.outfile = file.path(.outdir, 'genotype_pereira+er.tsv')

genotype_pereira %>>% write_df(.outfile)

#########1#########2#########3#########4#########5#########6#########7#########
## Normal samples from 1000 genomes

cachedir = '~/working/cancer'
kg_signif = file.path(cachedir, 'tidy.all.20130502.genotypes.topdrivers.signif.tsv.gz') %>>%
    read_tsv() %>>%
    dplyr::left_join(drivers, 'gene') %>>%
    dplyr::mutate(pathway= dplyr::coalesce(pathway, 'none')) %>>%
    (?.)

genotype_healthy = kg_signif %>>%
    dplyr::count(sample, pathway) %>>%
    tidyr::spread(pathway, n, 0L) %>>%
    dplyr::ungroup() %>>% dplyr::select(-sample) %>>%
    {ifelse(. > 0L, 1L, 0L)} %>>% as_tibble() %>>%
    dplyr::mutate(epsilon=1L) %>>%
    (?.)


#########1#########2#########3#########4#########5#########6#########7#########
## Make genotype matrix from GDC MAF

setwd('~/working/cancer')
.drivers = read_tsv('driver_genes.tsv') %>>% (?.)
top_drivers = .drivers %>>%
    dplyr::filter(refs >= 4) %>>%
    dplyr::rename(pathway=vogel_pathway) %>>%
    tidyr::replace_na(list(role='TSG', pathway='NA')) %>>%
    dplyr::select(gene, role, pathway) %>>%
    (?.)
init_genes = setNames(rep(1, nrow(top_drivers)), top_drivers$gene) %>>% (?.)

.pathway_genes = .drivers %>>%
    tidyr::replace_na(list(vogel_pathway='NA')) %>>%
    dplyr::mutate(pathway= str_split(vogel_pathway, ':')) %>>%
    tidyr::unnest() %>>%
    dplyr::count(pathway) %>>%
    (?.)
init_pathways = setNames(rep(1, nrow(.pathway_genes)), .pathway_genes$pathway) %>>% (?.)

maf_annot = read_tsv('gdc_topdrivers.maf.clean.tsv.gz') %>>%
    dplyr::left_join(top_drivers, by='gene') %>>%
    (?.)

maf_annot %>>%
    dplyr::count(cancer_type, sample) %>>%
    dplyr::count(cancer_type) %>>%
    dplyr::arrange(desc(nn))

maf_annot %>>%
    dplyr::mutate(width = end - start) %>>%
    ggplot(aes(width))+geom_histogram(center=0, binwidth=1)

#########1#########2#########3#########4

calc_mudens = function(.data) {.data %>>%
    tidyr::nest(-role, -pathway, -gene, -mutype) %>>%
    dplyr::select(-pathway) %>>%
    purrr::by_row(~{
        table(.$data[[1]]$start, dnn='pos') %>>%
        tibble::as_tibble() %>>%
        dplyr::mutate(pos=as.integer(pos),
          mut_density= n / sum(n), n=NULL)
    }) %>>%
    dplyr::select(-data) %>>%
    tidyr::unnest()
}

mut_density = maf_annot %>>%
    (? dplyr::count(., mutype)) %>>%
    dplyr::mutate(mutype= str_replace(mutype, 'splice', 'truncating')) %>>%
    dplyr::filter(mutype %in% c('missense', 'truncating')) %>>%
    calc_mudens() %>>% (?.)

mut_density %>>%
    dplyr::filter(mut_density > 0.1) %>>%
    dplyr::filter(!(role=='oncogene' & mutype=='truncating')) %>>%
    dplyr::arrange(role, gene, pos) %>>%
    (?.) %>>% wtl::write_df('recurrent_mutations.tsv.gz')

add_mudens_weight = function(.data) {
    .oncomut_density = mut_density %>>%
        dplyr::filter(role=='oncogene') %>>%
        dplyr::select(-role)
    left_join(.data, .oncomut_density, by=c('gene', 'pos', 'mutype')) %>>%
    dplyr::mutate(mudens_weight= ifelse(role=='oncogene', mut_density, 1)) %>>%
    tidyr::replace_na(list(mudens_weight=0)) %>>%
    dplyr::select(-mut_density)
}

#' @param tsg_missense numeric value relative to truncating mutation
summarise_mut_score_by_gene = function(.data, tsg_missense) {
    dplyr::mutate(.data, tsg_weight=
      ifelse(role!='oncogene' & mutype=='missense', tsg_missense, 1.0)) %>>%
    dplyr::group_by(sample, pathway, gene, role) %>>%
    dplyr::summarise(mut_score= sum(mudens_weight * tsg_weight)) %>>%
    dplyr::ungroup() %>>%
    dplyr::mutate(mut_score= pmin(1.0, mut_score))
}

maf_weighted = maf_annot %>>% dplyr::mutate(pos=start) %>>%
    add_mudens_weight() %>>%
    summarise_mut_score_by_gene(0.5) %>>%
    dplyr::arrange(sample, pathway, desc(mut_score)) %>>% (?.)
