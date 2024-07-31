#' @title reduceComps: Reduce the components of the model.
#' @param componentsOld The old components in the model.
#' @param pointsCopy Logical: is the spatial model for the points copy.
#' @param biasCopy Logical: is the bias model copy.
#' @param datasetName Name of the dataset of interest.
#' @param reducedTerms Terms to include in the new components.
#' @return New components.

reduceComps <- function(componentsOld, 
                        pointsCopy,
                        biasCopy,
                        datasetName,
                        reducedTerms) {
  
  all_comp_terms <- labels(terms(componentsOld))
  
  if (pointsCopy) {
    ##Make this function
    Main <- grepl('_spatial', all_comp_terms) & grepl('_field', all_comp_terms)
    Main <- sub("\\(.*", "", all_comp_terms[Main])
    
    if (paste0(datasetName,'_spatial') == Main) {
      
      Copy <-  grepl('_spatial', all_comp_terms) & grepl('copy', all_comp_terms) & gsub('\\(.*$', '', all_comp_terms) %in% reducedTerms
      
      compNext <- all_comp_terms[Copy][1]
      nextTerm <- sub("\\(.*", "", compNext)
      compNew <- paste0(nextTerm, '(main = geometry, model = ', datasetName, '_field)')
      all_comp_terms[Copy][1] <- compNew
      fixComps <- TRUE
      
    } else fixComps <- FALSE
    
  } else fixComps <- FALSE
  
  if (biasCopy) {
    
    MainBias <- grepl('_biasField', all_comp_terms) & grepl('_bias_field', all_comp_terms)
    MainBias <- sub("\\(.*", "", all_comp_terms[MainBias])
    
    if (paste0(datasetName,'_biasField') == MainBias) {
      
      CopyBias <-  grepl('_biasField', all_comp_terms) & grepl('copy', all_comp_terms)
      
      compBiasNext <- all_comp_terms[CopyBias][1]
      nextBiasTerm <- sub("\\(.*", "", compBiasNext)
      compBiasNew <- paste0(nextBiasTerm, '(main = geometry, model = ', datasetName, '_bias_field)')
      all_comp_terms[CopyBias][1] <- compBiasNew
      fixBiasComps <- TRUE
      
    } else fixBiasComps <- FALSE
    
  } else fixBiasComps <- FALSE
  
  comp_terms <- gsub('\\(.*$', '', all_comp_terms)
  
  comp_out <- comp_terms %in% reducedTerms
  
  if (all(comp_out)) reduced_components <- componentsOld
  else reduced_components <- update(componentsOld, paste0(' ~ . -', paste0(all_comp_terms[!comp_out], collapse = ' - ')))
  
  if (fixComps) {
    
    reduced_components <- update.formula(reduced_components, formula(paste0( '~ . -', compNext)))
    reduced_components <- update.formula(reduced_components, formula(paste0('~ . +', compNew)))
    
    if (sum(Copy) > 1) {
      newForm <- c()
      rmIndex <- which(Copy)[2:sum(Copy)]
      for (change in rmIndex) {
        
        newForm[change] <- gsub(paste0('copy = \"', Main, '\"'),
                                paste0('copy = \"', nextTerm, '\"'),
                                all_comp_terms[change])  
        
        
      }
      newForm <- na.omit(newForm)
      
      reduced_components <- update.formula(reduced_components, formula(paste0( '~ . -', paste0(all_comp_terms[rmIndex], collapse = '-'))))
      reduced_components <- update.formula(reduced_components, formula(paste0( '~ . +', paste0(newForm, collapse = '+'))))
      
    }
    
  }
  
  if (fixBiasComps) {
    
    reduced_components <- update.formula(reduced_components, formula(paste0( '~ . -', compBiasNext)))
    reduced_components <- update.formula(reduced_components, formula(paste0('~ . +', compBiasNew)))
    
    if (sum(CopyBias) > 1) {
      newForm <- c()
      rmIndex <- which(Copy)[2:sum(CopyBias)]
      for (change in rmIndex) {
        
        newForm[change] <- gsub(paste0('copy = \"', MainBias, '\"'),
                                paste0('copy = \"', nextBiasTerm, '\"'),
                                all_comp_terms[change])  
        
        
      }
      newForm <- na.omit(newForm)
      
      reduced_components <- update.formula(reduced_components, formula(paste0( '~ . -', paste0(all_comp_terms[rmIndex], collapse = '-'))))
      reduced_components <- update.formula(reduced_components, formula(paste0( '~ . +', paste0(newForm, collapse = '+'))))
      
    }
    
  }
  
  reduced_components
  
}