#' Outputs for bru_sdm_data
#' 

setClassUnion("listorNULL", c('list','NULL'))

setClass('bru_sdm_data', 
         slots = c(PO_data = 'listorNULL',
                   PA_data = 'listorNULL',
                   Mark_data = 'listorNULL',
                   ips = 'SpatialPointsDataFrame',
                   mesh = 'inla.mesh'))


setMethod('show', 'bru_sdm_data',
          function(x) {
            cat('\n')
            cat('Summary of inlabru_sdm data file:\n\n')
            
            if (!is.null(x@PA_data)) {
              
              cat('Present absence datasets:\n\n')
              PA_data <- data.frame(c('-----',names(x@PA_data)),
                                    c('',rep('|  ---  |', length(x@PA_data))),
                                    c('------------------',unlist(lapply(x@PA_data, function(dat) nrow(dat@coords)))))
              names(PA_data) <- c('Name:','','# of observations:')
              print.data.frame(PA_data[,1:3], row.names = FALSE, right = FALSE)
              cat('\n')
              
            } 
            #else {
            #    
            #  cat('Present absence datasets:\n\n')
            #  cat('No presenct absence datasets found.')
            #  cat('\n\n')
            #  
            #  }
            
            if (!is.null(x@PO_data)) {
              
              cat('Present only datasets:\n\n')
              PO_data <- data.frame(c('-----',names(x@PO_data)),
                                    c('',rep('|  ---  |', length(x@PO_data))),
                                    c('------------------',unlist(lapply(x@PO_data, function(dat) nrow(dat@coords)))))
              names(PO_data) <- c('Name:','','# of observations:')
              print.data.frame(PO_data[,1:3], row.names = FALSE, right = FALSE)
              cat('\n')
              
            }
            #else {
            #  
            #  cat('Present only datasets:\n\n')
            #  cat('No presenct only datasets found.')
            #  cat('\n\n')
            #  
            #  }
            
            
            
            if (!is.null(x@Mark_data)) {
              
              cat('Marks included:\n\n')
              mark_data <- data.frame(c('-----',unlist(lapply(x@Mark_data, function(dat) attributes(dat)$mark_name))),
                                      c('',rep('|  ---  |', length(x@Mark_data))),
                                      c('-----',unlist(lapply(x@Mark_data, function(dat) attributes(dat)$data_type))))
              #mark_data <- mark_data[!duplicated(lapply(x@Mark_data, function(dat) attributes(dat)$mark_name)),]
              mark_data <- unique(mark_data)
              names(mark_data) <- c('Name:','', 'Type:')
              print.data.frame(mark_data[,1:3], row.names = FALSE, right = FALSE)
              cat('\n')
              
              
            } 
            if (!is.null(attributes(x)$Timevariable)) {
              
              cat('Temporal variable included: ', attributes(x)$Timevariable)  
              
              
            }
            
          })
