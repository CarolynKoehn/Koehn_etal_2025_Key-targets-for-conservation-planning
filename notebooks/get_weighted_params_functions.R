#################################################
# Purpose: Functions for creating a table of weighted parameters for any data source based on NLCD classified parameters and NLCD overlap with data source classes
# Author: Carolyn Koehn
# Modified: 10/12/22
# Comments: This script creates a table of parameters to use in InVEST models.
#   I use weighted averages of NLCD-specific parameters that reflect the NLCD classes present in non-NLCD data classes.
#################################################

library(tidyverse)

create_comparison_table <- function(rasterStack, labels1, labels2) {
  #create data frame to store land class comparisons with layer 1 codes
  corrs <- data.frame(L1_code = terra::unique(rasterStack[[1]])[,1])
  #merge codes and class definitions
  corrs <- dplyr::left_join(corrs, labels1, by=c("L1_code" = "Code"))
  
  #nested loop to get number L2 pixels in each L1 class
  for (i in 1:nrow(corrs)) { #for each AFT class
    tab <- table(rasterStack[[2]][][rasterStack[[1]][]==corrs$L1_code[i]]) #create a table showing the distribution of NLCD classes in the AFT class
    for (j in 1:length(labels(tab)[[1]])) { #for the NLCD classes present in the AFT class
      for (code in 1:nrow(labels2)) { #compare to all NLCD classes
        if (labels(tab)[[1]][j] == as.character(labels2$Code[code])) { #match codes to find column name
          corrs[i,as.character(labels2$Code[code])] <- tab[[j]] #add to data frame how many pixels of NLCD are in the AFT class
        }
      }
    }
  }
  
  #add total area of each class
  corrs$Total_Area <- terra::freq(rasterStack[[1]])$count
  
  return(corrs)
}

#creates function that iterates through a row of the comparison table to get percent class present in each row
calculate_weights <- function(row.vec,code.col,label.col,totalArea.col,filterBelow) {
  total <- as.numeric(row.vec[totalArea.col])
  area_for_weights <- row.vec[-c(code.col,label.col,totalArea.col)]
  weights <- as.numeric(area_for_weights)/total
  names(weights) <- names(area_for_weights)
  return(weights)
}

get_weighted_params <- function(comparison.table, known.params,
                                code.col=1,label.col=2,
                                noData.col="No Data",
                                totalArea.col=ncol(comparison.table),
                                filterBelow=0.01) { #[fraction] if a land class is present at a level under the cutoff, it will not be included in the parameter calculation and its land area will be distributed to other land use classes
  #if any columns were given as names instead of indexes, convert here
  col.codes <- list(code.col=code.col,label.col=label.col,
                    noData.col=noData.col,totalArea.col=totalArea.col)
  col.codes <- lapply(col.codes,
                      FUN=function(x) {
                        if (!is.character(x)) {
                          x <- as.integer(x)
                        } else {
                          x <- which(colnames(comparison.table) == x)
                        }
                      })
  
  #create table to store final params
  params_table <- data.frame(lucode=comparison.table[,col.codes$code.col],
                             LULC_name=comparison.table[,col.codes$label.col])
  #add columns as specified by know params
  params_table[colnames(known.params)[3:ncol(known.params)]] <- NA
  
  #make comparison.table into matrix for apply()
  comp.table <- as.matrix(comparison.table)
  weights.array <- apply(comp.table, MARGIN=1,
                         FUN=function(x) calculate_weights(row.vec=x, code.col=col.codes$code.col, label.col=col.codes$label.col, totalArea.col=col.codes$totalArea.col,filterBelow=filterBelow))
  
  weights.mat <- as.matrix(weights.array)
  for(j in 1:ncol(weights.mat)) {
    if(weights.mat[noData.col,j]<filterBelow | is.na(weights.mat[noData.col,j])) {
      addUnder <- sum(weights.mat[weights.mat[,j]<filterBelow,j],na.rm=TRUE)
    } else {
      addUnder <- sum(weights.mat[weights.mat[,j]<filterBelow,j],na.rm=TRUE)
      addUnder <- addUnder+weights.mat[noData.col,j]
    }
    for(i in 1:nrow(weights.mat)) {
      if(weights.mat[i,j]<filterBelow | is.na(weights.mat[i,j]) | rownames(weights.mat)[i]==noData.col) {
        weights.mat[i,j] <- 0
      } else {
        weights.mat[i,j] <- weights.mat[i,j] + (weights.mat[i,j]*addUnder)
      }
    }
  }
  #since these weights don't quite add up to one, I'm adding the difference to the highest category. Shouldn't affect too much but makes the numbers add up
  for (j in 1:ncol(weights.mat)) {
    if(colSums(weights.mat)[j]!=1){
      weights.mat[which.max(weights.mat[,j]),j] <- weights.mat[which.max(weights.mat[,j]),j] + (1-colSums(weights.mat)[j])
    }
  }
  
  weights.table <- data.frame(lucode = comparison.table[,col.codes$code.col],
                              LULC_name = comparison.table[,col.codes$label.col])
  weights.table <- cbind(weights.table,t(weights.mat))
  
  #getting final weighted params from weights.table and known.params
  #i want to do it by param
  for(p in 3:ncol(params_table)) {
    wparam_table <- weights.table
    for(i in 1:nrow(known.params)) {
      wparam_table[,as.character(known.params[i,1])] <- wparam_table[,as.character(known.params[i,1])]*known.params[i,p]
      wparam_vec <- apply(wparam_table[,3:ncol(wparam_table)],MARGIN=1,FUN=sum)
    }
    params_table[,p] <- wparam_vec
  }
  
  return(params_table)
}