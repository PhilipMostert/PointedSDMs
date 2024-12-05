#' @description Function to make the datasets into likelihoods.
#' @title \emph{makeLhoods}: function to make likelihoods.
#' @param data A list of sf objects containing the datasets for which likelihoods need to be constructed.
#' @param formula A list of formulas to add to the likelihoods.
#' @param family A list of vectors containing the families within each dataset.
#' @param mesh An \code{fm_mesh_2d} object.
#' @param ips Integration points used.
#' @param paresp The response variable name for the presence absence datasets.
#' @param ntrialsvar The trials variable name for the presence absence datasets.
#' @param markstrialsvar The trial variable name for the binomial marks.
#' @param speciesname The name of the species variable used.
#' @param speciesindex A vector containing the numeric index of where the species occurs in the data
#' @param samplers A list of integration domains for the datasets.
#' @param pointcovs A vector of the point covariates used in the model.

makeLhoods <- function(data, formula, family, mesh, ips,
                       paresp, ntrialsvar, markstrialsvar,
                       speciesname, speciesindex, samplers,
                       pointcovs = NULL) {
  
  Likelihoods <- list()

  for (dataset in 1:length(data)) {
    
    for (species in 1:length(data[[dataset]])) {
      
      if (length(data[[dataset]][[species]]) != 0) {
      
      Ntrialsvar <- list()
      
      if (!is.null(ntrialsvar)) {
        
        if (ntrialsvar %in% names(data[[dataset]][[species]])) {
          
          Ntrialsvar[[1]] <- data.frame(data[[dataset]][[species]])[,ntrialsvar]
          
        } else Ntrialsvar[[1]] <- NA
        
      } else Ntrialsvar[[1]] <- NA
      
      if (!is.null(markstrialsvar)) {
        
        if (markstrialsvar %in% names(data[[dataset]][[species]])) {
          
          Ntrialsvar[[2]] <- data.frame(data[[dataset]][[species]])[,markstrialsvar]
          
        } else Ntrialsvar[[2]] <- NA
        
      } else Ntrialsvar[[2]] <- NA
      
      for (process in 1:length(family[[dataset]])) {
        
        Likindex <- length(Likelihoods) + 1
        
        if (family[[dataset]][process] == 'binomial') {
         
          if (!all(is.na(Ntrialsvar[[1]])) || !all(is.na(Ntrialsvar[[2]]))) {
            
            if (as.character(formula[[dataset]][[species]][[process]][['LHS']])[2] == paresp) Ntrials <- Ntrialsvar[[1]]
            
            else 
              if (!any(is.na(Ntrialsvar[[2]]))) Ntrials <- Ntrialsvar[[2]]
              else Ntrials <- 1
            
          } else Ntrials <- 1
        } 
        else Ntrials <- 1
        
        if (!is.null(speciesname)) {
          
          if (length(speciesindex) != 0) {
            
            speciesRep <- data.frame(rep(unique(data.frame(data[[dataset]][[species]])[,speciesname]), nrow(ips)))
            names(speciesRep) <- speciesname
            speciesRep$speciesSpatialGroup <- speciesRep[,speciesname]
            IPS <- ips
            
            namesKeep <- names(IPS)[names(IPS) %in% c('weight', '.block',names(data[[dataset]][[species]]))]
            IPS <- IPS[, namesKeep]
            IPS <- cbind(IPS, speciesRep)
            
          }
          
        }
        else 
          if (!is.null(samplers[[names(data)[[dataset]]]])) {
            
            IPS <- NULL
            
          }
        else IPS <- ips
        
        if (!is.null(pointcovs)) {
        
        if (family[[dataset]][process] == 'cp' && any(pointcovs %in% names(data[[dataset]][[species]]))) {
          
          #pointcovsIn <- pointcovs[pointcovs %in%  names(data[[dataset]][[species]])]
          #formula[[dataset]][[species]][[process]][['LHS']] <- reformulate(deparse(formula[[dataset]][[species]][[process]][['LHS']][[3]]), paste0('geometry + ', pointcovsIn))
        
        }
        }

        # bru_like_list will use the like-tag for list names; inlabru >= 2.12.0
        if (is.null(names(data[[dataset]])[species])) nameGive <- names(data)[[dataset]]
        else nameGive <- names(data[[dataset]])[species]
        
        like_name <- paste0(nameGive, '_', sub(' .*', '', as.character(formula[[dataset]][[species]][[process]][['LHS']])[2]))
        
        Likelihoods[[Likindex]] <- inlabru::like(formula = formula[[dataset]][[species]][[process]][['LHS']], ## but obs change these in function call
                                                 include = formula[[dataset]][[species]][[process]][['RHS']],
                                                 data = data[[dataset]][[species]], 
                                                 Ntrials = Ntrials,
                                                 ips = IPS,
                                                 domain = list(geometry = mesh),
                                                 samplers = samplers[[names(data)[[dataset]]]],
                                                 family = family[[dataset]][process],
                                                 tag = like_name)
        
        names(Likelihoods)[[Likindex]] <- like_name
        
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