library(tidyverse)
loadNamespace('cowplot')

.labels0 = 'z,x,y,fill
AA,3.5,3.5,#888888
AB,8,3.5,#e41a1c
AC,10,3.5,#377eb8
BA,3.5,8,#e41a1c
BB,8,8,#888888
BC,10,8,#4daf4a
CA,3.5,10,#377eb8
CB,8,10,#4daf4a
CC,10,10,#888888
' %>% read_csv()

.tbl0 = tidyr::crossing(x=seq_len(10L), y=x) %>%
    dplyr::mutate(
      first = ifelse(y <= 6L, 'A', ifelse(y <= 9L, 'B', 'C')),
      second= ifelse(x <= 6L, 'A', ifelse(x <= 9L, 'B', 'C')),
      z= paste0(first, second)) %>%
    dplyr::left_join(.labels0 %>% dplyr::select(z, fill), by='z') %>%
    print()

.p0 = ggplot(.tbl0, aes(x, y))+
    geom_tile(aes(fill=fill), show.legend=FALSE)+
    scale_fill_identity()+
    geom_hline(data=tibble(y=c(6.5, 9.5)), aes(yintercept=y))+
    geom_vline(data=tibble(x=c(6.5, 9.5)), aes(xintercept=x))+
    geom_text(data=.labels0, aes(label=z), size=6, colour='#ffffff')+
    scale_y_reverse()+
    coord_fixed()+
    theme_void()
.p0

.labels1 = 'z,x,y,fill
AB,16,12,#e41a1c
AC,36,12,#377eb8
BA,18,30.5,#e41a1c
BC,37.5,30.5,#4daf4a
CA,13.5,38.5,#377eb8
CB,33.5,38.5,#4daf4a
' %>% read_csv()

.tbl1 = tidyr::crossing(x=seq_len(40L), y=x) %>%
    dplyr::mutate(
      first = ifelse(y <= 24L, 'A', ifelse(y <= 36L, 'B', 'C')),
      second = ifelse(first == 'A', ifelse(x <= 30L, 'B', 'C'),
               ifelse(first == 'B', ifelse(x <= 40 * 6 / 7, 'A', 'C'),
                                    ifelse(x <= 40 * 2 / 3, 'A', 'B'))),
      z= paste0(first, second)) %>%
    dplyr::left_join(.labels1 %>% dplyr::select(z, fill), by='z') %>%
    print()

.p1 = ggplot(.tbl1, aes(x, y))+
    geom_tile(aes(fill=fill), show.legend=FALSE)+
    scale_fill_identity()+
    geom_hline(data=tibble(y=c(24.5, 36.5)), aes(yintercept=y))+
    geom_text(data=.labels1, aes(label=z), size=6, colour='#ffffff')+
    scale_y_reverse()+
    coord_fixed()+
    theme_void()
.p1

.plt = cowplot::plot_grid(.p0, .p1, nrow=1, labels=c('Current', 'Proposed'))
.plt
ggsave('mutation_order_prob.png', .plt, width=8, height=4.2)

# #######1#########2#########3#########4#########5#########6#########7#########

p = c(0.5, 0.3, 0.2)

sum(2 * p[1] * p[2] / (1 - sum(p ** 2)),
    2 * p[2] * p[3] / (1 - sum(p ** 2)),
    2 * p[3] * p[1] / (1 - sum(p ** 2)))

sum(p[1] * p[2] / (1 - p[1]) + p[2] * p[1] / (1 - p[2]),
    p[2] * p[3] / (1 - p[2]) + p[3] * p[2] / (1 - p[3]),
    p[3] * p[1] / (1 - p[3]) + p[1] * p[3] / (1 - p[1]))

purrr::rerun(3L, sample(seq_along(p), 2, replace=FALSE, prob=p))
