#---------------------------------------------------------------------#
# Project: Investigating the contribution of residual unexplaind 
#          variability components in nonlinear mixed effect approach 
# Program: Generate simulated data using mrgsolve
# Author: Mutaz M. Jaber <jaber038@umn.edu> 
# Date created: 9/5/21
# Date modified: 9/17/21
#---------------------------------------------------------------------#
library(mrgsolve) 

args  <- commandArgs(trailingOnly=TRUE) 
nsubs <- as.numeric(args[1])
nsim  <- as.numeric(args[2])
TDOSE <- as.numeric(args[3])
DOSE  <- as.numeric(args[4])
TYPE  <- args[5] 
BASE  <- args[6]
PER   <- args[7] 

MODEL <- mrgsolve::mread(BASE) 



if (PER == 'A1') {
        sigmat <- diag(c(0.0025, 1e-04))
} else if (PER == 'A2') {
        sigmat <- diag(c(0.0025, 4e-4))
} else if (PER == 'A3' | PER == 'All') {
        sigmat <- diag(c(0.0025, 0.0016))
} else {
        sigmat <- diag(c(0.0025, 0.0))
} 
  
model <- mrgsolve::smat(MODEL, sigmat)

  
GetData <- function(model, 
                    nsim=100,
                    TDOSE,
                    DOSE,
                    SAMPLE,
                    nsubs,
                    TYPE,
                    per=c('B','A1', 'A2', 'A3', 'S1','SL1', 'SL2', 'SL3', 'S2', 'D','TD1','TD2','All')
                    ) {

        if (inherits(per, "character")) {
                per <- match.arg(per) 
        }

        if (inherits(design, "character")) {
                design <- match.arg(design) 
        } 
        
	sims <- list()
	dose <- list()
	conc <- list()
	base <- list()
	
        BNAME <- c("ID", "TIME", "DV", "AMT", "MDV", "EVID") # Not sure if CMT should be added at this point
        if (TYPE == 'int') {
                SAMPLE <- c(0.5, 1, 2, 4, 6, 8, 12, 24, 48) 
        } else if (TYPE == 'spa') {
                SAMPLE <- c(2, 12, 24, 48) 
        } else if (TYPE == 'spa1') {
                SAMPLE <- c(2, 6, 12, 24, 48)
        } else if (TYPE == 'spa2') {
                SAMPLE <- s(1, 2, 6, 12, 24, 48)
        }
        } else {
                stop('Study design must be declared (int or spa)')
        }

        if (per =='B' | per == 'A1' | per == 'A2' | per == 'A3') {
        	event <- ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$time %in% SAMPLE,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                
                                colnames(base[[i]]) <- BNAME
                              }

      } else if (per == 'S1') {
               event <- expand.ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                 		set.seed(i)                  
                 		SAM <- c(0, SAMPLE)
                 		dl <- purrr::map(event$ID, ~ sample(SAMPLE + rnorm(length(SAMPLE), 0, 5/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID) 
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			                    dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				                dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                colnames(base[[i]]) <- BNAME
                                base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
	       }
      } else if (per == 'SL1') {
               event <- expand.ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                 		set.seed(i) 
                 		SAM <- c(0, SAMPLE)
                 		dl <- purrr::map(event$ID, ~ sample(SAMPLE + rnorm(length(SAMPLE), 0, 10/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID) 
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				    dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                colnames(base[[i]]) <- BNAME
                                base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
                                base[[i]]$TIME <- ifelse(base[[i]]$AMT >0 & base[[i]]$TIME !=0,
                                                         0, base[[i]]$TIME
                                                         )
                                base[[i]] <- base[[i]][!(base[[i]]$DV ==0 & base[[i]]$EVID==0),]
	       }
      } else if (per == 'SL2') {
               event <- expand.ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                 		set.seed(i) 
                 		SAM <- c(0, SAMPLE)
                 		dl <- purrr::map(event$ID, ~ sample(SAMPLE + rnorm(length(SAMPLE), 0, 15/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID) 
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                colnames(base[[i]]) <- BNAME
                                base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
                                base[[i]]$TIME <- ifelse(base[[i]]$AMT >0 & base[[i]]$TIME !=0,
                                                         0, base[[i]]$TIME
                                                         )
                                base[[i]] <- base[[i]][!(base[[i]]$DV ==0 & base[[i]]$EVID==0),]
	       }
      } else if (per == 'SL3') {
               event <- expand.ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                 		set.seed(i) 
                 		SAM <- c(0, SAMPLE)
                 		dl <- purrr::map(event$ID, ~ sample(SAMPLE + rnorm(length(SAMPLE), 0, 30/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID) 
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                colnames(base[[i]]) <- BNAME
                                base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
                                base[[i]]$TIME <- ifelse(base[[i]]$AMT >0 & base[[i]]$TIME !=0,
                                                         0, base[[i]]$TIME
                                                         )
                                base[[i]] <- base[[i]][!(base[[i]]$DV ==0 & base[[i]]$EVID==0),]
	       }
	} else if (per == 'S2') {
	       event <- expand.ev(amt=DOSE, time=TDOSE, ID=1:nsubs)
               for (i in 1:nsim) {
                                set.seed(i) 
                 		SAM <- c(0, SAMPLE)
                 		dl <- purrr::map(event$ID, ~ sample(SAMPLE + runif(length(SAMPLE), -5/60, 5/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID)
                 		set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                colnames(base[[i]]) <- BNAME
                                base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
	       }
	} else if (per == 'TD1') {
	      
	                      for (i in 1:nsim) {
	                        set.seed(i) 
	                        TDOSE = rnorm(nsubs, 0, 5/60)
	       			event <- expand.ev(amt=DOSE, time=TDOSE) 
	       			
	       			set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt"))
			        dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				dose[[i]]$Cc <- 0 
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$time %in% SAMPLE,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                
                                colnames(base[[i]]) <- BNAME
                                base[[i]][base[[i]]$AMT != 0, ]$TIME <- 0 
                              }
        } else if (per =='TD2') {

	                      for (i in 1:nsim) {
	                        set.seed(i) 
	                        TDOSE = rnorm(nsubs, 0, 10/60)
	       			        event <- expand.ev(amt=DOSE, time=TDOSE) 
	       			        set.seed(i) 
                            sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt"))
			                dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
				            dose[[i]]$Cc <- 0 
                                # Dose data 
                            conc[[i]] <- sims[[i]][sims[[i]]$time %in% SAMPLE,]
                                # Arrange data

                            base[[i]] <- rbind(dose[[i]], conc[[i]])
                            base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                
                            colnames(base[[i]]) <- BNAME
                            base[[i]][base[[i]]$AMT != 0, ]$TIME <- 0 
                              }             
           
        } else if (per == 'D') {
        	
               for (i in 1:nsim) {
                 		 set.seed(i) 
		                 DOSE = rnorm(nsubs, DOSE, 12)
		                 event <- ev(amt=DOSE, time=TDOSE, ID=1:nsubs, replicate=FALSE)
		                 
		                 set.seed(i) 
                                sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt"))
                                dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
                                dose[[i]]$Cc <- 0 
                                dose[[i]]$amt <- 120
                                # Dose data 
                                conc[[i]] <- sims[[i]][sims[[i]]$time %in% SAMPLE,]
                                # Arrange data

                                base[[i]] <- rbind(dose[[i]], conc[[i]])
                                base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                                
                                colnames(base[[i]]) <- BNAME
                              }
               
        } else if (per == 'All') {
          	
               for (i in 1:nsim) {
                 		 
                 		set.seed(i); TDOSE = rnorm(nsubs, 0, 5/60)
         			    set.seed(i); DOSE = rnorm(nsubs, DOSE, 12)
          			    event <- ev(amt=DOSE, time=TDOSE, ID=1:nsubs, replicate=FALSE)
                 		SAM <- c(0, SAMPLE)
                 		set.seed(i); dl <- purrr::map(event$ID, ~ sample(SAMPLE + rnorm(length(SAMPLE), 0, 5/60), length(SAMPLE)))
                 		idata <- dplyr::select(event, ID) 
				        set.seed(i) 
                        sims[[i]] <- as.data.frame(mrgsim(model, event, outvars="Cc", carry_out="amt,evid,cmt", idata=idata, deslist=dl, descol="ID"))
			            dose[[i]] <- sims[[i]][sims[[i]]$amt != 0,]
                        dose[[i]]$Cc <- 0 
                        dose[[i]]$amt <- 120
                                # Dose data 
                        conc[[i]] <- sims[[i]][sims[[i]]$amt==0,]
                                # Arrange data

                        base[[i]] <- rbind(dose[[i]], conc[[i]])
                        base[[i]] <- subset(base[[i]][order(base[[i]]$ID, base[[i]]$time),],select=c(ID,time,Cc,amt,cmt,evid))
                        colnames(base[[i]]) <- BNAME
                        base[[i]]$TIME <- rep(SAM, length=length(unique(base[[i]]$ID)))
	       }          
        }
               
	names(base) <- paste0("BASE", seq_along(base))


        if (!dir.exists(paste("../../data", TYPE, per, sep="/"))) {
                dir.create(paste("../../data",TYPE,per,sep="/"), recursive=TRUE)
        }
	for (i in 1:nsim) {
	    NAM <- paste("../../data",TYPE,PER,paste0("dat",i,".csv"),sep='/')
	    write.csv(base[[i]], NAM,quote=FALSE,row.names=F,na='.')
	 }

}

GetData(model,nsim,TDOSE,DOSE,SAMPLE,nsubs,TYPE,PER) 
