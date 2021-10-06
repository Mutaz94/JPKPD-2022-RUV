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

# # Generating plots 
plot.1 <- ggplot(DAT[DAT$X=='RUV',], aes(per, log(Mean), group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        geom_errorbar(aes(ymin=log(X95CIlo), ymax=log(X95CIup)), width=0.2)+
        labs(color='Study design', x='Perturbation', y =expression(log(Deviation)), title=expression(RUV))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'M', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(M), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 

plot.2 <- ggplot(DAT[DAT$X=='RUV',], aes(per, rBias*100, group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        labs(color='Study design', x='Perturbation', y =expression(rBias), title=expression(RUV))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 
plot.3 <- ggplot(DAT[DAT$X=='RUV',], aes(per, rRMSE*100, group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        labs(color='Study design', x='Perturbation', y =expression(rRMSE), title=expression(RUV))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(labels=c(expression(SD[8]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 

plot.4 <- ggplot(DAT[DAT$X=='CL',], aes(per, log(Mean), group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        geom_errorbar(aes(ymin=log(X95CIlo), ymax=log(X95CIup)), width=0.2)+
        labs(color='Study design', x='Perturbation', y =expression(log(Deviation)), title=expression(CL))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 
plot.5 <- ggplot(DAT[DAT$X=='V',], aes(per, log(Mean), group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        geom_errorbar(aes(ymin=log(X95CIlo), ymax=log(X95CIup)), width=0.2)+
        labs(color='Study design', x='Perturbation', y =expression(log(Deviation)), title=expression(V))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 
plot.6 <- ggplot(DAT[DAT$X=='Q',], aes(per, log(Mean), group=design))+
        geom_point(aes(color=design), size=1)+
        geom_line(aes(color=design), size=.75)+
        geom_errorbar(aes(ymin=log(X95CIlo), ymax=log(X95CIup)), width=0.2)+
        labs(color='Study design', x='Perturbation', y =expression(log(Deviation)), title=expression(Q))+
        scale_x_discrete(limits=c('B', 'A1', 'A2', 'A3', 'D', 'S1', 'SL1', 'SL2', 'SL3','TD1', 'TD2', 'All'), 
        labels=c('Base', expression(A[1]), expression(A[2]), expression(A[3]), expression(D[v]), 
                  expression(S[5]), expression(S[10]), expression(S[15]), expression(S[30]), 
                  expression(D[t5]), expression(D[t10]), expression(All)))+
        scale_color_discrete(limits=c('SD1', 'SD3', 'SD2', 'SD4'),
                             labels=c(expression(SD[9]), expression(SD[6]), expression(SD[5]), expression(SD[4]))) 
pdf('RUV-JPKPD-Figures.pdf', width=12, height=8)
plot.1
plot.2
plot.3
plot.4
plot.5
plot.6
dev.off()
