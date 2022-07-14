#' @description Function to make the datasets into likelihoods.
#' @title \emph{makeLhoods}: function to make likelihoods.
#' @param data A list of sp objects containing the datasets for which likelihoods need to be constructed.
#' @param formula A list of formulas to add to the likelihoods.
#' @param family A list of vectors containing the families within each dataset.
#' @param mesh An inla.mesh object.
#' @param ips Integration points used.
#' @param paresp The response variable name for the presence absence datasets.
#' @param ntrialsvar The trials variable name for the presence absence datasets.
#' @param markstrialsvar The trial variable name for the binomial marks.
#' @param speciesname The name of the species variable used.
#' @param speciesindex A vector containing the numeric index of where the species occurs in the data
#' @param samplers A list of integration domains for the datasets.

makeLhoods <- function(data, formula, family, mesh, ips,
                       paresp, ntrialsvar, markstrialsvar,
                       speciesname, speciesindex, samplers) {
  
  Likelihoods <- list()

  for (dataset in 1:length(data)) {
    
    for (species in 1:length(data[[dataset]])) {
      
      if (length(data[[dataset]][[species]]) != 0) {
      
      Ntrialsvar <- list()
      
      if (!is.null(ntrialsvar)) {
        
        if (ntrialsvar %in% names(data[[dataset]][[species]]@data)) {
          
          Ntrialsvar[[1]] <- data[[dataset]][[species]]@data[,ntrialsvar]
          
        } else Ntrialsvar[[1]] <- NA
        
      } else Ntrialsvar[[1]] <- NA
      
      if (!is.null(markstrialsvar)) {
        
        if (markstrialsvar %in% names(data[[dataset]][[species]]@data)) {
          
          Ntrialsvar[[2]] <- data[[dataset]][[species]]@data[,markstrialsvar]
          
        } else Ntrialsvar[[2]] <- NA
        
      } else Ntrialsvar[[2]] <- NA
      
      for (process in 1:length(family[[dataset]])) {
        
        Likindex <- length(Likelihoods) + 1
        
        if (family[[dataset]][process] == 'binomial') {
         
          if (!all(is.na(Ntrialsvar[[1]])) || !all(is.na(Ntrialsvar[[2]]))) {
            
            if (as.character(formula[[dataset]][[species]][[process]][['LHS']])[2] == paresp) Ntrials <- Ntrialsvar[[1]]
            
            else Ntrials <- Ntrialsvar[[2]]
            
          } else Ntrials <- 1
        } 
        else Ntrials <- 1
        
        if (!is.null(speciesname)) {
          
          if (length(speciesindex) != 0) {
            
            speciesRep <- data.frame(rep(unique(data[[dataset]][[species]]@data[,speciesname]), nrow(ips@coords)))
            names(speciesRep) <- speciesname
            
            IPS <- ips
            IPS@data <- cbind(ips@data, speciesRep)
            
          }
          
        }
        else IPS <- ips
        
        Likelihoods[[Likindex]] <- inlabru::like(formula = formula[[dataset]][[species]][[process]][['LHS']], ## but obs change these in function call
                                                 include = formula[[dataset]][[species]][[process]][['RHS']],
                                                 data = data[[dataset]][[species]], 
                                                 Ntrials = Ntrials,
                                                 mesh = mesh,
                                                 ips = IPS,
                                                 samplers = samplers[[names(data)[[dataset]]]],
                                                 family = family[[dataset]][process])
        
        if (is.null(names(data[[dataset]])[species])) nameGive <- names(data)[[dataset]]
        else nameGive <- names(data[[dataset]])[species]
        
        names(Likelihoods)[[Likindex]] <- paste0(nameGive, '_', as.character(formula[[dataset]][[species]][[process]][['LHS']])[2])
        
      } 
      
      } else {
       
        Likindex <- length(Likelihoods) + 1
         
        Likelihoods[[Likindex]] <- NULL
        
      }
        
    }
    
  }
  
  Likelihoods <- Likelihoods[!unlist(lapply(Likelihoods, is.null))]
  
  Likelihoods  
  
}