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

.outdir = '~/Dropbox/working/tumor'
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
