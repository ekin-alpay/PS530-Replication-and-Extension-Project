############################################################################
## Determinants of Sanctions Effectiveness
## Replication file - Functions
## Authors: Navin Bapat, Tobias Heinrich, Yoshiharu Kobayashi, Cliff Morgan, 
## Date: September 8, 2012
############################################################################



# Function that runs all combinations of the variable for one imputed dataset
# Inputs are the x-th imputed dataset from the list of imputed datasets (complete.datasets); comb is the matrix of all combinations of
# variables; the.DV specifies the name of the dependent variable used in the estimations; the.vars is a vector of strings that give all
# possible predictors to be included. 
# x denotes the number of the imputed datasets
run.mods <- function(x, complete.datasets, comb, the.DV, cases., the.vars, security)
{  
  data <- complete.datasets[[x]]
  
  if(security == 1) data <- subset(data, security == 1)
  
  Y <- data[,which(colnames(data) == the.DV)]  # Extract the DV
  X <- data[,which(colnames(data) %in% the.vars)] ## Extracts all predictors to be used
  var.names <- colnames(X)
  
  run.single.model <- function(picked.vars, cases.)
  { ## Function that takes one of the combinations and estimates the model,
    ## and stores coefficients, SEs, and log likelihood
    if(cases. == "imposition" | cases. == "threatonly")
    {                  
      picked.cases <- which(complete.datasets[[x]][,which(colnames(complete.datasets[[x]]) == cases.)]==1)
      X.temp <- as.matrix(X[picked.cases,picked.vars])
      Y.temp <- Y[picked.cases]
    }
    
    if(cases. == "NA")
    {
      X.temp <- as.matrix(X[,picked.vars])
      Y.temp <- Y
    }
    
    mod <- glm(Y.temp ~ 1 + X.temp, family=binomial(link="probit"))
    
    out <- matrix(NA, 2, length(the.vars))
    
    out[1,picked.vars] <- coef(mod)[-1]
    out[2,picked.vars] <- sqrt(diag(vcov(mod))[-1])
    return(out)
  }
  
  estimated <- aaply(.data=comb, .margins=2, .fun=run.single.model, cases.=cases., .progress="none")  # Evaluates run.single.model() for each of the combinations
  if(sum(comb[,1]) == sum(seq(1:length(the.vars))))
  {
    estimated2 <- array(NA, dim=c(1,2, length(the.vars)))
    estimated2[1,,] <- estimated
    estimated <- estimated2
  }
  dimnames(estimated)[[2]] <- c("Coefficients", "SEs")
  dimnames(estimated)[[3]] <- var.names
  return(estimated)
}


## This function takes a list of analysis output from different imputed data sets 
## and brings that output together
openandclean <- function(filee, the.vars)
{
  load(filee)
  
  tmp <- strsplit(x=filee, split="Included")[[1]][2]
  n.vars <- as.numeric(gsub("\\D", "", tmp))
  
  n.imp <- length(results)
  get.out <- function(x, what)
  {
    tmp <- melt(results[[x]][,what,])
    if(ncol(tmp) == 1) tmp <- cbind(rownames(tmp), tmp)
    if(ncol(tmp) == 3) tmp <- tmp[which(is.na(tmp[,3])==FALSE),][,c(2,3)]
    tmp <- data.frame(tmp, what)
    return(tmp)    
  }
  
  all.coefs <- adply(.data=seq(1:n.imp), .margins=1, .fun=get.out, what="Coefficients")
  all.SEs   <- adply(.data=seq(1:n.imp), .margins=1, .fun=get.out, what="SEs")
  all.T <- all.coefs
  all.T$value <- abs(as.matrix(all.coefs$value)/as.matrix(all.SEs$value))
  all.T$what <- "t Statistics (Absolute)"
  
  temp <- rbind(all.coefs, all.T)
  colnames(temp) <- c("Imputation", "Variable", "Value", "Type")
  temp <- merge(temp, proper.names, by="Variable")
  temp$Imputation <- NULL
  temp$Variable <- NULL
  temp$Nvars <- n.vars
  return(temp)
}
