.indir = '~/Dropbox/working/cancer/genotypes'
.infile = file.path(.indir, 'genotype_pereira+er.tsv')
genotype_pereira = read_tsv(.infile) %>>% (?.)

.names = names(genotype_pereira)
.freqs = genotype_pereira %>>% colMeans() %>>% (. / sum(.))

.nested = genotype_pereira %>>%
    dplyr::group_by(s=rowSums(.)) %>>%
    tidyr::nest() %>>%
    dplyr::arrange(s) %>>%
    (?.)

.nested$data[[3]]
