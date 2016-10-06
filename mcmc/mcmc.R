library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(logical=FALSE))

#########1#########2#########3#########4#########5#########6#########7#########

.model = rstan::stan_model('simplex.stan')

working_dir = '~/Dropbox/working/cancer/genotypes'
.infiles = list.files(working_dir, 'tsv$', full.names=TRUE)

make_data = function(.infile, threshold = 0.3, intercept = 0.1, sd_coef = 0.2) {
    .genotypes = readr::read_tsv(.infile)
    list(
    threshold = threshold,
    intercept = intercept,
    sd_coef = sd_coef,
    n_samples = nrow(.genotypes),
    n_coefs = ncol(.genotypes),
    genotypes = .genotypes)
}
.data = make_data(.infiles[1])

.fit = rstan::sampling(.model, data=.data, iter=1000, chains=3)

#########1#########2#########3#########4

print(.fit)
summary(.fit)

stan_plot(.fit, pars=c('coefs'))+
    coord_cartesian(xlim=c(0, 1))

stan_trace(.fit, pars=c('coefs'))+
    coord_cartesian(ylim=c(0, 1))

stan_hist(.fit, pars=c('coefs'))+
    coord_cartesian(xlim=c(0, 1))

stan_dens(.fit, pars=c('coefs'), color=NA, separate_chains=TRUE)+
    coord_cartesian(xlim=c(0, 1))

#########1#########2#########3#########4#########5#########6#########7#########

.results = wtl::map_par(.infiles, ~{
    .base = tools::file_path_sans_ext(basename(.x))
    .outfile = paste0(.base, '.stanfit.rds')
    message(.outfile)
    .fit = rstan::sampling(.model, make_data(.x), iter=10000, chains=3)
    saveRDS(.fit, .outfile)
    .fit
}, .cores=2, .cluster='PSOCK')
