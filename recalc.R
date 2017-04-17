.annot = list(
    pathway= c('A', 'B'),
    annotation= c('0011', '1100'))

.s1 = .annot
.s1$sample = c('0001', '0010', '0100', '1000')
.s1

.s2 = .annot
.s2$sample = c('0011', '0101', '1001', '0110', '1010', '1100')
.s2

.s3 = .annot
.s3$sample = c('0111', '1011', '1101', '1110')
.s3

.tablize = function(.data) {
    tibble(genotype= .data[['sample']]) %>>%
      dplyr::mutate(
        g_A= str_extract(genotype, '\\d\\d$'),
        g_B= str_extract(genotype, '^\\d\\d'),
        g_A1= str_extract(g_A, '\\d$') %>>% as.integer(),
        g_A2= str_extract(g_A, '^\\d') %>>% as.integer(),
        g_B1= str_extract(g_B, '\\d$') %>>% as.integer(),
        g_B2= str_extract(g_B, '^\\d') %>>% as.integer(),
        s_A= str_count(g_A, '1'),
        s_B= str_count(g_B, '1'),
        a_A= pmax(s_A - 1L, 0L),
        a_B= pmax(s_B - 1L, 0L))
}
.tablize(.s1)
.tablize(.s2)
.tablize(.s3)

.calc_p_raw = function(.data, .w_gene, .th_path) {
    .data %>>% dplyr::mutate(
        p_base= (.w_gene[['g_A1']] ^ g_A1) * (.w_gene[['g_A2']] ^ g_A2) * (.w_gene[['g_B1']] ^ g_B1) * (.w_gene[['g_B2']] ^ g_B2),
        excl= (.th_path[['A']] ^ a_A) * (.th_path[['B']] ^ a_B),
        p_raw= p_base * excl)
}
.w_gene = seq_len(4) %>>% (. /sum(.)) %>>% setNames(sprintf('g_%s%d', rep(c('A', 'B'), each=2), rep(c(1, 2), 2)))
.calc_p_raw(.tablize(.s1), .w_gene, .th_path)
.calc_p_raw(.tablize(.s2), .w_gene, .th_path)
.calc_p_raw(.tablize(.s3), .w_gene, .th_path)

.cal_loglik = function(.data, .th_path) {
    .tbl = .tablize(.data)
    .w_gene = .tbl %>>% dplyr::summarise_at(vars(matches('g_\\w\\d')), sum) %>>% (. / sum(.)) %>>% (?.)
    .D = .tablize(.s3) %>>% .calc_p_raw(.w_gene, .th_path) %>>% (p_raw) %>>% sum()
    .tbl %>>%
        .calc_p_raw(.w_gene, .th_path) %>>%
        dplyr::mutate(p= p_raw / .D, loglik= log(p))
}

.testcase = .annot
.testcase$sample = c('0011', '0101', '1001', '0110', '1010', '1100', '1100')
.testcase

.th_path = c(A=1.2, B=1.2)
.cal_loglik(.testcase, .th_path) %>>% (?.) %>>% (sum(.$loglik))
.cal_loglik(.s3, .th_path) %>>% (?.) %>>% (sum(.$loglik))

system2('likeligrid', c('-s4'), input=jsonlite::toJSON(.testcase))

.con = gzfile('test.json.gz', 'w')
jsonlite::write_json(.testcase, .con, pretty=TRUE)
close(.con)
