

# Walter Reinisch, Vivek Pradhan, Saira Ahmad, Zhen Zhang, Jeremy D Gale, Alternative endoscopy reading paradigms determine score reliability and effect size in ulcerative colitis, Journal of Crohn's and Colitis, 2023;, jjad134, https://doi.org/10.1093/ecco-jcc/jjad134

rm(list=ls())
require(haven)
fnam <- 'mcs_all.RData'
reRun <- FALSE

setwd(dirWork <- './')
sasDir <- "C:\\Program Files\\SASHome\\SASFoundation\\9.4"
sasCmd <- paste0('sas.exe -PRINT ', dirWork, '\\cdrsb.lst -LOG ',
                 dirWork,'\\cdrsb.log -SYSIN ', dirWork, '\\runModel.sas')

dat <- as.data.frame(read_sas('mcs_all.sas7bdat'))
mnam <- c('LR','Super-reader','CR1','CR2','ConCR','ConLCR1','ConLCR2')
n <- length(unique(dat$METHOD))
nmodel <- length(modNam <- paste0(rep(c('Normal','OrdLogit'),each=2), rep(c('',' w/ covariate'), 2)))

if(!file.exists(fnam) || reRun){
  out <- oc <- numeric();  rnam <- character()
  ptm <- proc.time()[3] 
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      mat <- droplevels(dat[dat$METHODN %in% c(i, j),]) 
      mat <- mat[order(mat$WINDOW), ]
      write.csv(file='foo.csv', mat[,c('WINDOW','SUBJID','MCSFLXN','ANTITNFN')], row.names=FALSE, na='')
      
      for(k in 1:nmodel) if(file.exists(fnam1 <- paste0('icc',k,'.csv'))) file.remove(fnam1)
      for(k in 1:nmodel) if(file.exists(fnam1 <- paste0('fit',k,'.csv'))) file.remove(fnam1)
      setwd(sasDir);  return_code <- shell(sasCmd);  setwd(dirWork) 
      #if(return_code!=0) stop('model running issues.')

      rnam <- c(rnam, paste(mnam[i], 'vs.', mnam[j]))
      for(k in 1:nmodel){
        icc <- read.csv(paste0('icc',k,'.csv'))
        vec <- format(round( as.matrix(icc[match(icc$WINDOW, c(0,12)), 
                                           c('Estimate','Lower','Upper')]),  digits=3 ), nsmall=3)
        out <- rbind(out, data.frame(
                       'Baseline'=paste0(vec[1,1], ' (', vec[1,2], ',', vec[1,3] ,')'), 
                       'Week 12'=paste0(vec[2,1], ' (', vec[2,2], ',', vec[2,3] ,')')
               ))
        
        fit <- read.csv(paste0('fit',k,'.csv'))
        tmp <- as.data.frame(t(fit$Value))
        names(tmp) <- gsub(' \\(smaller is better\\)', '', paste(fit$WINDOW, fit$Descr, sep=':'))
        oc <- rbind(oc, tmp)
      }
      
    }
  }
  cputime <- as.numeric(proc.time()[3]-ptm)/60 
  cat('\nCPUtime', round(cputime,3), 'minutes: completed!','\n')
  
  rnamFull <- paste0( rep(rnam, each=nmodel), ': ', rep(modNam, length(rnam)) )
  row.names(out) <- rnamFull;   row.names(oc) <- rnamFull
  names(out) <- gsub('\\.', ' ', names(out))
  print(out)
  save(file=fnam, out, oc)
}else load(fnam)

ind <- function(modelName) return(seq(which(modNam == modelName), nrow(oc), by=nmodel))
dif <- oc[ind('Normal'),] - oc[ind('OrdLogit'),]
# dif <- oc[ind('Normal w/ covariate'),] - oc[ind('OrdLogit w/ covariate'),]
cat(paste0(100*round(mean(dif[,c('0:AIC','12:AIC')] > 0),4), 
           '% chance that ordinal logit model has a smaller AIC (better) than Normal model.\n')) 

dif <- oc[ind('OrdLogit'),] - oc[ind('OrdLogit w/ covariate'),]
cat(paste0(100*round(mean(dif[,c('0:AIC','12:AIC')] > 0),4), 
           '% chance that ordinal logit model w/ covariates has a smaller AIC (better) than w/o covariates.\n')) 


#=============== backup
# # R version
library(lmerTest)
require(reshape2)
m <- lmer(MCSFLXN ~ (1|SUBJID), REML=FALSE, data=mat)
#Create function to calculate ICC from fitted model
calc.icc <- function(m) {
  sumy <- summary(m)
  (sumy$varcor$SUBJID[1]) / (sumy$varcor$SUBJID[1] + sumy$sigma^2)
}
calc.icc(m)
set.seed(123)
boot.icc <- bootMer(m, calc.icc, nsim=1000)
#Draw from the bootstrap distribution the usual 95% upper and lower confidence limits
quantile(boot.icc$t, c(0.025, 0.975))

# not run
















