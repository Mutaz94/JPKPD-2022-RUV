#---------------------------------------------------------------------#
# Project: Investigating the contribution of residual unexplaind 
#          variability components in nonlinear mixed effect approach 
# Program: Generate simulated data 
# Author: Mutaz M. Jaber <jaber038@umn.edu> 
# Date created: 9/5/21
# Date modified: 9/5/21
#---------------------------------------------------------------------#

args  <- commandArgs(tailingOnly=TRUE) 
nsubs <- args[1]
nsim  <- args[2]
TDOSE <- args[3]
DOES  <- args[4]
TYPE  <- args[5] 
LEVEL <- args[6] 
BASE  <- args[7]
PER   <- args[8] 

MODEL <- mrgsolve::mread(BASE) 

if (TYPE == 'int') {
        SAMPLE <- c(0.5, 1, 2, 4, 6, 8, 12, 24, 48) 
} else if (TYPE == 'spa') {
        SAMPLE <- c(2, 12, 24, 48) 
} else {
        stop('Study design must be declared (int or spa)')
}

if (LEVEL == 1) {
        sigmat <- diag(c(0.0025, 1e-04))
} else if (LEVEL == 2) {
        sigmat <- diag(c(0.0025, 4e-4))
} else if (LEVEL == 3) {
        sigmat <- diag(c(0.0025, 0.0016))
} else {
        sigmat <- diag(c(0.0025, 0.0))
} 

  
model <- mrgsolve::smat(MODEL, sigmat)

  
GetData <- function(model, 
                    nsim=100,
                    per=c('B','A1', 'A2', 'A3', 'S1', 'D', 'TD','All'),
                    design=c('INT', 'SPA')) {

        if (inherits(per, "character")) {
                per <- match.arg(per) 
        }

        if (inherits(design, "character")) {
                design <- match.arg(design) 
        } 

        BNAME <- c("ID", "TIME", "DV", "AMT", "MDV", "EVID") # Not sure if CMT should be added at this point

        if (per == 'B') {
               event <- ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars='Cc', carry_out='amt,evid,cmt'))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                
                                colnames(base[[i]]) <- BNAME
		    }
	}
	names(base) <- paste0("BASE", seq_along(base))
	for (i in 1:nsim) {
	    NAM[i] <- paste("../../data",type,PER,'dat',i,".csv",sep='/')
	    write.csv(base[[i]], NAM[i],quote=FALSE,row.names=F,na='.')
	 }










