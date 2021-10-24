#-----------------------------------------------------------------
# Project: Investigating the contribution of residual unexplained 
#          variability components in NLME approach
# Author: Mutaz M. Jaber
# Date created: 9.25.21
# Date modified: 9.25.21
#-----------------------------------------------------------------

library(ggplot2) 
library(ggthemes)
theme_set(theme_base())



files <- list.files('../../results/', pattern='csv') 
dat <- list()

for (j in files) {
        dat[[j]] <- read.csv(paste0('../../results/',j))
}

names(dat) <- gsub("*.csv", "", names(dat)) 

DAT <- purrr::map_df(dat, ~as.data.frame(.x), .id='id')
DAT <- tidyr::separate(DAT,id, c('design', 'per'), "-")
ID <- unique(DAT$id)

PK_PARAM_LIST = unique(DAT$X)
# # # Generating plots 

for (PK in PK_PARAM_LIST) {
    assign(paste0('plot.', PK), ggplot(DAT[DAT$X==PK,], aes(per, Mean, group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        geom_errorbar(aes(ymin=X95CIlo, ymax=X95CIup), width=0.2)+
        labs(color='Study design', x='Perturbation', y =expression(log(Deviation)), title=PK)+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3','D', 'D2', 'M', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]),expression(D[v2]), expression(M),
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4])))) 


    assign(paste0('rBias.', PK), ggplot(DAT[DAT$X==PK,], aes(per, rBias*100, group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        labs(color='Study design', x='Perturbation', y =expression(rBias), title=PK)+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3','D', 'D2', 'M', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]),expression(D[v2]), expression(M),
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4])))) 

    assign(paste0('rRMSE.', PK), ggplot(DAT[DAT$X==PK,], aes(per, rRMSE*100, group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        labs(color='Study design', x='Perturbation', y =expression(rRMSE), title=PK)+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3','D', 'D2', 'M', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]),expression(D[v2]), expression(M),
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4])))) 

}


 pdf('Dev-JPKPD-Figures.pdf', width=12, height=8)
        mget(ls(pattern='plot'))
 dev.off()

 pdf('rBias-JPKPD-Figures.pdf', width=12, height=8)
        mget(ls(pattern='rBias'))
 dev.off()

 pdf('rRMSE-JPKPD-Figures.pdf', width=12, height=8)
        mget(ls(pattern='rRMSE'))
 dev.off()
