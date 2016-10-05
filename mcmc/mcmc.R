library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#########1#########2#########3#########4#########5#########6#########7#########

working_dir = '~/Dropbox/working/tumor'
.infile = file.path(working_dir, 'genotype_pereira+er.tsv')
genotype_pereira = readr::read_tsv(.infile)

.data = list(
    threshold = 0.3,
    intercept = 0.1,
    sd_coef = 0.2,
    cancer_genotypes = genotype_pereira)
.data$n_coefs= ncol(.data$cancer_genotypes)
.data$n_samples= nrow(.data$cancer_genotypes)

.model = rstan::stan_model('simplex.stan')
.fit = rstan::sampling(.model, data=.data, iter=1000, chains=3)

print(.fit)
summary(.fit)
plot(.fit, pars=c('coefs'))+coord_cartesian(xlim=c(0, 1))
# pairs(.fit)
rstan::traceplot(.fit)

rstan::stan_trace(.fit, pars=c('coefs'))+
    coord_cartesian(ylim=c(0, 1))

rstan::stan_hist(.fit, pars=c('coefs'))+
    coord_cartesian(xlim=c(0, 1))

rstan::stan_dens(.fit, pars=c('coefs'), separate_chains=TRUE)+
    coord_cartesian(xlim=c(0, 1))
