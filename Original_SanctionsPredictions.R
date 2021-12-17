############################################################################
## Determinants of Sanctions Effectiveness
## Replication file 
## Authors: Navin Bapat, Tobias Heinrich, Yoshiharu Kobayashi, Cliff Morgan, 
## Date: September 8, 2012
############################################################################

rm(list=ls())

set.seed(1981)

## The following are options for the analysis
## Set the working directory; default is the current directory.
## This is where the 
wkdir <- "/Users/th5/Dropbox/Projects/SanxPredV2/estimation/"

## If reimpute=TRUE, then multiple imputation is run anew.
## If reimpute=FALSE, then existing files will be loaded from
## the filename "CompleteDatasets.Rdata". If they do
## not exist, the replication will not be executed.
reimpute <- TRUE

## This sets whether multiple cores of your computer should be used.
## doParallel=TRUE is unlikely to work under Windows as the multicore R
## package is not natively supported. It will work under OS X Lion or 
## later. Modify the parallel options below (such as number of cores, etc.)
## Thus: If you use Windows, please set doParallel=FALSE. 
doParallel <- TRUE

## Each graph in the analysis is based on about 10 million data points. 
## As this may cause severe memory issues for some computers, you take a random
## draw of the estimates. (In trials, we found that even just 10% of the 
## simulation provided a very accurate approximation to distribution of 
## all simulations.)
## randomFrac sets the fraction of simulations to be randomly drawn.
randomFrac <- 1


#####################
## Getting started ##
#####################
## Loading packages
library(Amelia)
library(foreign)
library(plyr)
library(ggplot2)
library(reshape)


## Set directories
setwd(wkdir)
dir.create(paste(getwd(), "/saved-runs/", sep=""), showWarnings=FALSE)
dir.create(paste(getwd(), "/figures/", sep=""), showWarnings=FALSE)


## Parallel options | doParallel=TRUE
if(doParallel == TRUE)
{
  library(foreach)
  library(doMC)
  library(multicore)
  registerDoMC(cores=3)
}


## Loading data ##
data.raw <- read.dta(file="SancDataForAnalysisVer7.dta") 

data.raw <- get(load("/Users/apple/Desktop/Replication-1/SancDataForAnalysisVer7.RData"))
data.raw


##data.raw <- data.raw[-which(colnames(data.raw) %in% c("issue11",  "ti1", 
#                                                    "sender1name", "targetname", "endyear", "sender1", "targetstate","sancstartyear"))]
##This command was given as a command in the syntax, as the dataset loads only with the required numerical variables. I did not exclude this anyway. 


# Define cases that were only threats
data.raw$threatonly <- 1 - data.raw$imposition

# Define security-related sanctions
data.raw$security <- 0
data.raw$security[data.raw$issue_trade == 0 & data.raw$issue_econ == 0 & data.raw$issue_env == 0] <- 1

# Drop unused variables
data.raw <- data.raw[, c("ally", "capratio", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue",
                         "multisender", "rivalry", "asendercost", "smartsanc", "CI", "salience",
                         "Tdependence", "financial", "US", "csendercost", "ctargetcost", "atargetcost", "sendercost",
                         "targetcost", "LRsuccess", "Tsuccess", "threatonly", "imposition",
                         "Slnpop", "Tlnpop", "lntrade", "Rsuccess", "issue_trade", "issue_env",
                         "issue_econ", "TRsuccess", "security")]



# Transforming data
to.be.binary <- c("SOECD", "TOECD", "Sinterstate", "Sinternal", "Tinterstate", "Tinternal")
for(i in 1:length(to.be.binary))
{
  k <- which(colnames(data.raw) == to.be.binary[i])
  data.raw[,k] <- as.integer(data.raw[,k])-1
}


## Standardize several variables
to.be.stand <- c("CI", "Tlnpop", "Slnpop", "lntrade", "capratio")
for(i in 1:length(to.be.stand))
{
  k <- which(colnames(data.raw) == to.be.stand[i])
  data.raw[,k] <- (data.raw[,k] - mean(data.raw[,k], na.rm=T))/sd(data.raw[,k], na.rm=T)
}    


# Proper Names
proper.names <- rbind(c("lntrade", "Trade"), c("capratio", "Capability Ratio"), c("rivalry", "Rivalry"),
                      c("ally", "Ally"), c("asendercost", "Sender Costs"), c("atargetcost", "Target Costs"),
                      c("institution", "IO Involvement"), c("multisender", "Multiple Senders"), c("multiissue", "Multiple Issues"),
                      c("import", "Import Restrictions"), c("carrot", "Carrots"), c("carrotimp", "Carrots"), c("smartsanc", "Smart Sanctions"),
                      c("Tdem", "Democratic Target"), c("Sdem", "Democratic Sender"), c("CI", "Target Instability"),
                      c("executive", "Threat by Executive"), c("bspecific", "Behavior Specificity"), c("commit", "Sender Commitment"),
                      c("targetcost", "Target Costs"), c("sendercost", "Sender Costs"), c("ctargetcost", "Target Costs"),
                      c("csendercost", "Sender Costs"), c("export","Export Restrictions"), c("financial", "Financial Sanctions"),
                      c("US", "United States Sender"), c("salience", "High Issue Saliency"), c("Tdependence", "Target Trade Dependence"),
                      c("threatonly", "Threat-Only Cases"), c("imposition", "Imposed-Only Cases"), c("lntrade", "Trade Volume"),
                      c("Tlnpop", "Target Population"), c("Slnpop", "Sender Population"), c("LRsuccess", "Success Outcome"), c("Tsuccess", "Success Outcome"))
colnames(proper.names) <- c("Variable", "VariableProper")
proper.names <- as.data.frame(proper.names)


## Missingness graphs

library(dplyr)
library(plyr)
the.vars <- c("ally", "capratio", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue",
              "multisender", "rivalry", "sendercost", "smartsanc", "CI", "targetcost", "Tdependence",
              "salience", "financial", "US", "ctargetcost", "csendercost", "atargetcost", 
              "asendercost", "LRsuccess", "Tsuccess")
data.raw$Groups <- "All"
data.raw$Groups[data.raw$threatonly==1] <- "Threat Only"
data.raw$Groups[data.raw$imposition==1] <- "Imposed Only"
data.ms1 <- ddply(.data=data.raw, .variables="Groups", 
                  .fun=function(y) adply(.data=y, .margins=2, .fun=function(x) sum(is.na(x))/nrow(y)))
data.ms2 <- adply(.data=data.raw, .margins=2, .fun=function(x) sum(is.na(x))/nrow(data.raw))
data.ms <- rbind(data.ms1, data.frame(Groups="All", X1=data.ms2[,1], V1=data.ms2[,2]))
data.ms <- subset(data.ms, X1 %in% the.vars)
colnames(data.ms) <- c("Groups", "Variable", "Missingness")
data.ms <- merge(data.ms, proper.names, by="Variable")
data.ms <- data.ms[(data.ms$Variable == "Tsuccess" & data.ms$Groups == "All") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "Tsuccess" & data.ms$Groups == "Imposed Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "LRsuccesxss" & data.ms$Groups == "Threat Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "targetcost" & data.ms$Groups == "Threat Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "targetcost" & data.ms$Groups == "All") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "sendercost" & data.ms$Groups == "Threat Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "sendercost" & data.ms$Groups == "All") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "ctargetcost" & data.ms$Groups == "Threat Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "ctargetcost" & data.ms$Groups == "Imposed Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "csendercost" & data.ms$Groups == "Threat Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "csendercost" & data.ms$Groups == "Imposed Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "atargetcost" & data.ms$Groups == "Imposed Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "atargetcost" & data.ms$Groups == "All") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "asendercost" & data.ms$Groups == "Imposed Only") == FALSE,]
data.ms <- data.ms[(data.ms$Variable == "asendercost" & data.ms$Groups == "All") == FALSE,]
data.ms <- subset(data.ms, Variable != "threatonly")
data.ms <- subset(data.ms, Variable != "imposition")
data.ms$VariableProper <- as.character(data.ms$VariableProper)

## Graph missingness
g <- ggplot(data=data.ms, aes(x=VariableProper, y=Missingness))
g <- g + geom_point(size=4) + coord_flip() + geom_hline(yintercept=c(0,1), size=1)
g <- g + opts(panel.grid.major = theme_line(size=.051), panel.grid.minor = theme_blank(), panel.background=theme_rect(fill="grey90", colour="grey90"), axis.line=theme_segment(colour=NA), strip.text.y = theme_text())
g <- g + facet_wrap(~ Groups, nrow=1) + scale_x_discrete("") + scale_y_continuous("Percentage Missing")
ggsave(file=paste(getwd(), "/Figures/Missingness.pdf", sep=""), width=12, height=11)
rm(data.ms,g)


# Do Multiple Imputations
install.packages("mlr")
library("mlr")
install.packages("Amelia")
library("Amelia")
if(reimpute == TRUE)
{
  imputations <- amelia(data.raw, m=6,
                        idvars=c("threatonly", "imposition", "Groups"),
                        ords=c("ally", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue", "multisender",  
                               "rivalry", "asendercost", "smartsanc", "salience", "financial", "US", "csendercost", "ctargetcost",
                               "atargetcost", "sendercost", "targetcost", "LRsuccess", "Tsuccess", "Rsuccess", "TRsuccess",
                               "security"))
  write.amelia(obj=imputations, file.stem = "impdata_test", format="dta")
  complete.datasets <- imputations$imputations
  save(complete.datasets, file="CompleteDatasets.Rdata")
}
if(reimpute == FALSE) load("CompleteDatasets.Rdata")



######################
## DEFINE FUNCTIONS ##
######################
## The sourced script loads several functions used for the sensitivity analysis
source("SanctionsPredictions Functions.R")



################################
## THE ACTUAL RUNS BEGIN HERE ##
################################
## Define the runs
the.DVs <- c("Restrictive", "Less Restrictive")
models <- list(One=list(name="All Cases",
                        the.vars <- c("ally", "capratio", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue",
                                      "multisender", "rivalry", "csendercost", "smartsanc", "CI", "ctargetcost", "Tdependence",
                                      "salience", "financial", "US")),
               Two=list(name="Threats Only",
                        the.vars <- c("ally", "capratio", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue",
                                      "multisender", "rivalry", "asendercost", "smartsanc", "CI", "atargetcost", "Tdependence",
                                      "salience", "financial", "US")),
               Three=list(name="Impositions Only",
                          the.vars <- c("ally", "capratio", "carrot", "Sdem", "Tdem", "export", "institution", "multiissue",
                                        "multisender", "rivalry", "sendercost", "smartsanc", "CI", "targetcost", "Tdependence",
                                        "salience", "financial", "US")))
samples <- c("Security", "All")                          


tm <- proc.time()[3]
count <- 0
## Loop that runs all estimations
for(j in 1:length(the.DVs))
{
  for(m in 1:length(samples))
  {
    for(k in 1:length(models))
    {
      the.vars <- models[[k]][[2]]
      if(models[[k]]$name == "Impositions Only") cases <- "imposition"
      if(models[[k]]$name == "All Cases") cases <- "NA"
      if(models[[k]]$name == "Threats Only") cases <- "threatonly"
      
      if(samples[m] == "Security") security <- 1
      if(samples[m] == "All") security <- 0
      
      if(the.DVs[j] == "Restrictive")
      {
        if(cases == "imposition" | cases == "NA") the.DV <- "Rsuccess"
        if(cases == "threatonly") the.DV <- "TRsuccess"
      }        
      
      if(the.DVs[j] == "Less Restrictive")
      {
        if(cases == "imposition" | cases == "NA") the.DV <- "LRsuccess"
        if(cases == "threatonly") the.DV <- "Tsuccess"
      }        
      
      dir.create(paste(getwd(), "/saved-runs/", sep=""), showWarnings=FALSE)
      dir.name <- paste(models[[k]]$name, "-", samples[m], "-", the.DVs[j], "/", sep="")
      dir.create(paste(getwd(), "/saved-runs/", dir.name, sep=""), showWarnings=FALSE)        
      
      ## The actual estimations
      cat("\nDone for this specification (out of 18): ")
      for(i in 1:length(the.vars))
      {
        tmp <- file.exists(paste(getwd(), "/saved-runs/", dir.name, "Included", i, ".Rdata", sep=""))
        comb <- combn(x=length(the.vars), m=i)
        if(tmp == FALSE)
        {
          results <- alply(.data=seq(1:length(complete.datasets)), .margins=1,
                           .fun=run.mods, complete.datasets=complete.datasets, comb=comb, 
                           the.DV=the.DV, cases.=cases, the.vars=the.vars, security=security,
                           .parallel=doParallel, .progress="none")
          save(results, file=paste(getwd(), "/saved-runs/", dir.name, "Included", i, ".Rdata", sep=""))
          rm(results)
        }
        cat(i," ", sep="")
      }
      cat("\n")
      
      # Combine the estimates
      tmp <- file.exists(paste(getwd(), "/saved-runs/", dir.name, "AllCombined.Rdata", sep=""))
      cat("Combining the estimates for graphing.\n")
      if(tmp == FALSE)
      {
        file.list <- c()
        for(i in 1:length(the.vars)) file.list <- c(file.list, paste(getwd(), "/saved-runs/", dir.name, "Included", i, ".Rdata", sep=""))
        out <- ldply(.data=file.list, .fun=openandclean, .parallel=doParallel, the.vars=the.vars, .progress="none")
        save(out, file=paste(getwd(), "/saved-runs/", dir.name, "AllCombined.Rdata", sep=""))
        rm(out)
      }
      dir.name <- paste(models[[k]]$name, "-", samples[m], "-", the.DVs[j], "/", sep="")
      load(file=paste(getwd(), "/saved-runs/", dir.name, "AllCombined.Rdata", sep=""))
      
      ## Plot Densities of coefficients/ t stats after combining the estimates
      graph.name <- paste(models[[k]]$name, "-", samples[m], "-", the.DVs[j], sep="")
      ind <- sample(1:nrow(out), size=round(randomFrac*nrow(out)), replace=FALSE)
      g <- ggplot(data=out[ind,], aes(x=Value, y=..scaled..))
      g <- g + geom_density(size=.15, fill="grey20", alpha=1/3)
      g <- g + facet_grid(VariableProper ~ Type, scales="free")
      g <- g + opts(panel.grid.major = theme_blank(), panel.grid.minor = theme_blank(), panel.background=theme_rect(fill="grey90", colour="grey90"), axis.line=theme_segment(colour=NA), strip.text.y = theme_text())
      g <- g + scale_y_continuous(breaks=NULL, labels=NA, limits=c(0,1.05)) + labs(x="", y="")
      g <- g + geom_line(data=data.frame(Type="Coefficients", x=c(0,0), y=c(0,1.05), Imputation=1), aes(x=x, y=y), size=.6, colour="gray40")
      g <- g + geom_line(data=data.frame(Type="t Statistics (Absolute)", x=c(1.65,1.65), y=c(0,1.05), Imputation=1), aes(x=x, y=y), size=.6, colour="gray40")
      g <- g + geom_point(data=data.frame(Type="t Statistics (Absolute)", x=c(0,9), y=c(0,0), Imputation=1), aes(x=x, y=y), size=.8, colour=NA)
      g <- g + geom_point(data=data.frame(Type="Coefficients", x=c(-.6,1.5), y=c(0,0), Imputation=1), aes(x=x, y=y), size=.8, colour=NA)
        
      graph.name <- paste(models[[k]]$name, "-", samples[m], "-", the.DVs[j], sep="")
      ggsave(file=paste(getwd(), "//saved-runs/", dir.name, "/", graph.name, ".pdf", sep=""), width=12, height=17)
      rm(out, g, ind)
      cat(paste("Finishing ", graph.name, "\n", sep=""))
      
      ## Time 
      cat("\n \n")
      thus.far <- proc.time()[3] - tm
      cat(round(thus.far/60),"Minutes have elapsed.")
      count <- count + 1
      to.go <- thus.far/count*12
      cat(" Roughly",round(to.go, di=0), "Minutes to go.")
      cat("\n \n")
    }
  }
}



####################################################################################
## Copying the figures that were generated into a separate folder, renaming them  ##
## to what they are named in the actual paper and in the web appendix
####################################################################################
dir.create(paste(getwd(), "/figures/", sep=""), showWarnings=FALSE)

# Figures for the paper
file.copy(from=paste(getwd(), "/saved-runs/All Cases-All-Less Restrictive/All Cases-All-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/Figure1.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Impositions Only-All-Less Restrictive/Impositions Only-All-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/Figure2.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Threats Only-All-Less Restrictive/Threats Only-All-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/Figure3.pdf", sep=""))

# Figures for the web appendix
file.copy(from=paste(getwd(), "/saved-runs/All Cases-Security-Less Restrictive/All Cases-Security-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA1.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Impositions Only-Security-Less Restrictive/Impositions Only-Security-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA2.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Threats Only-Security-Less Restrictive/Threats Only-Security-Less Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA3.pdf", sep=""))

file.copy(from=paste(getwd(), "/saved-runs/All Cases-All-Restrictive/All Cases-All-Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA4.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Impositions Only-All-Restrictive/Impositions Only-All-Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA5.pdf", sep=""))
file.copy(from=paste(getwd(), "/saved-runs/Threats Only-All-Restrictive/Threats Only-All-Restrictive.pdf", sep=""),
          to=paste(getwd(), "/figures/FigureA6.pdf", sep=""))

cat("\n \n \n \n \n")
cat("You have now replicated 'Determinants of Sanctions Effectiveness' \nby Bapat, Heinrich, Kobayashi, and Morgan.\n")
cat("The figures that correspond to those in the paper and the web \nappendix can be found in the 'figures' directory.")
cat("\n \n \n \n \n")
