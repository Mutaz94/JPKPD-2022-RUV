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



files <- list.files('../../results/') 
dat <- list()
for (j in files) {
        dat[[j]] <- read.csv(paste0('../../results/',j))
}

names(dat) <- gsub("*.csv", "", names(dat)) 

DAT <- purrr::map_df(dat, ~as.data.frame(.x), .id='id')
DAT <- tidyr::separate(DAT,id, c('design', 'per'), "-")
ID <- unique(DAT$id)

# Generating plots 
plot.1 <- ggplot(DAT[DAT$X=='RUV',], aes(per, log(Mean), group=design))+
        geom_point(aes(color=design))+
        geom_line()+
        geom_errorbar(aes(ymin=log(X95CIlo), ymax=log(X95CIup)))+
        labs(x='Perturbation', y =expression(log(Deviation)), title=expression(RUV))
pdf('RUV.pdf', height=10, width=12)
print(plot.1)
dev.off()
