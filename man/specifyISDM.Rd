% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specifyISDM.R
\name{specifyISDM}
\alias{specifyISDM}
\title{R6 class for creating a \code{startISDM} object.}
\description{
A data object containing the data and the relevant information about the integrated model. The function \code{\link{startISDM}} acts as a wrapper in creating one of these objects. The output of this object has additional functions within the object which allow for further specification and customization of the integrated model.
}
\examples{

## ------------------------------------------------
## Method `specifyISDM$plot`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 library(ggplot2)
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                            Projection = proj, 
                            responsePA = 'Present')
 
  #Create plot of data
  organizedData$plot()

}
}

## ------------------------------------------------
## Method `specifyISDM$addBias`
## ------------------------------------------------

 \dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                             Projection = proj, 
                             responsePA = 'Present')
 
#Add bias field to eBird records
organizedData$addBias(datasetNames = 'eBird')

}
}

## ------------------------------------------------
## Method `specifyISDM$updateFormula`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, 
                           responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Remove Forest from eBird
 organizedData$updateFormula(datasetName = 'eBird', Formula = ~ . - Forest)
 
 #Add some scaling to Forest for Parks
 organizedData$updateFormula('Parks', newFormula = ~ I(. +(Forest+1e-6)*scaling))
 
 #Now dd scaling to components
 organizedData$changeComponents(addComponent = 'scaling') 
 
}
}

## ------------------------------------------------
## Method `specifyISDM$changeComponents`
## ------------------------------------------------

\dontrun{

 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh, 
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Remove Forest from components
 organizedData$changeComponents(removeComponent = 'Forest')

}

}

## ------------------------------------------------
## Method `specifyISDM$priorsFixed`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Add prior to Forest
 organizedData$priorsFixed(Effect = 'Forest', mean.linear = 2, prec.linear = 0.1)

}
}

## ------------------------------------------------
## Method `specifyISDM$specifySpatial`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Specify the shared spatial field
 organizedData$specifySpatial(sharedSpatial = TRUE,
                       prior.range = c(1,0.001),
                       prior.sigma = c(1,0.001))

} 
}

## ------------------------------------------------
## Method `specifyISDM$changeLink`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh, 
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Specify the shared spatial field
 organizedData$changeLink('Parks', 'logit')
 
 
} 
}

## ------------------------------------------------
## Method `specifyISDM$spatialBlock`
## ------------------------------------------------

 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Specify the spatial block
 organizedData$spatialBlock(k = 2, rows = 2, cols = 1, plot = FALSE)
 
} 


## ------------------------------------------------
## Method `specifyISDM$addSamplers`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           Projection = proj, 
                           responsePA = 'Present')
 
#Add integration domain for the eBird records
organizedData$addSamplers(datasetName = 'eBird', Samplers = SolitaryTinamou$region)

}
}

## ------------------------------------------------
## Method `specifyISDM$specifyRandom`
## ------------------------------------------------

\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           Projection = proj, 
                           responsePA = 'Present',
                           pointsSpatial = copy)
 
#Add integration domain for the eBird records
organizedData$specifyRandom(copyModel =  list(beta = list(fixed = TRUE)))

}
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-specifyISDM-help}{\code{specifyISDM$help()}}
\item \href{#method-specifyISDM-print}{\code{specifyISDM$print()}}
\item \href{#method-specifyISDM-plot}{\code{specifyISDM$plot()}}
\item \href{#method-specifyISDM-addBias}{\code{specifyISDM$addBias()}}
\item \href{#method-specifyISDM-updateFormula}{\code{specifyISDM$updateFormula()}}
\item \href{#method-specifyISDM-changeComponents}{\code{specifyISDM$changeComponents()}}
\item \href{#method-specifyISDM-priorsFixed}{\code{specifyISDM$priorsFixed()}}
\item \href{#method-specifyISDM-specifySpatial}{\code{specifyISDM$specifySpatial()}}
\item \href{#method-specifyISDM-changeLink}{\code{specifyISDM$changeLink()}}
\item \href{#method-specifyISDM-spatialBlock}{\code{specifyISDM$spatialBlock()}}
\item \href{#method-specifyISDM-addSamplers}{\code{specifyISDM$addSamplers()}}
\item \href{#method-specifyISDM-specifyRandom}{\code{specifyISDM$specifyRandom()}}
\item \href{#method-specifyISDM-new}{\code{specifyISDM$new()}}
\item \href{#method-specifyISDM-samplingBias}{\code{specifyISDM$samplingBias()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-help"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-help}{}}}
\subsection{Method \code{help()}}{
Function to provide documentation for a \code{specifyISDM} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$help(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Not used}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Documentation.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-print"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-print}{}}}
\subsection{Method \code{print()}}{
Prints the datasets, their data type and the number of observations, as well as the marks and their respective families.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$print(...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{...}}{Not used.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-plot"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-plot}{}}}
\subsection{Method \code{plot()}}{
Makes a plot of the points surrounded by the boundary of the region where they were collected. The points may either be plotted based on which dataset they come from, or which species group they are part of (if \code{speciesName} is non-\code{NULL} in \code{\link{intModel}}).
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$plot(datasetNames, Boundary = TRUE, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetNames}}{Name of the datasets to plot. If this argument is missing, the function will plot all the data available to the model.}

\item{\code{Boundary}}{Logical: should a boundary (created using the \code{Mesh} object) be used in the plot. Defaults to \code{TRUE}.}

\item{\code{...}}{Not used.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A ggplot object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 library(ggplot2)
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                            Projection = proj, 
                            responsePA = 'Present')
 
  #Create plot of data
  organizedData$plot()

}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-addBias"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-addBias}{}}}
\subsection{Method \code{addBias()}}{
Function used to add additional spatial fields (called \emph{bias fields}) to a selected dataset present in the integrated model. \emph{Bias fields} are typically used to account for sampling biases in opportunistic citizen science data in the absence of any covariate to do such.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$addBias(
  datasetNames = NULL,
  allPO = FALSE,
  biasField = NULL,
  copyModel = TRUE,
  shareModel = FALSE,
  temporalModel = list(model = "ar1")
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetNames}}{A vector of dataset names (class \code{character}) for which a bias field needs to be added to. If \code{NULL} (default), then \code{allPO} has to be \code{TRUE}.}

\item{\code{allPO}}{Logical: should a bias field be added to all datasets classified as presence only in the integrated model. Defaults to \code{FALSE}.}

\item{\code{biasField}}{An \code{inla.spde} object used to describe the bias field. Defaults to \code{NULL} which uses \code{\link[INLA]{inla.spde2.matern}} to create a Matern model for the field.}

\item{\code{copyModel}}{Create copy models for all the of the datasets specified with either \code{datasetNames} or \code{allPO}. The first dataset in the vector will have its own spatial effect, and the other datasets will "copy" the effect with shared hyperparameters. Defaults to \code{TRUE}.}

\item{\code{shareModel}}{Share a bias field across the datasets specified with \code{datasetNames}. Defaults to \code{FALSE}.}

\item{\code{temporalModel}}{List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details. \code{temporalName} needs to be specified in \code{intModel} prior.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A bias field to the model.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ \dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                             Projection = proj, 
                             responsePA = 'Present')
 
#Add bias field to eBird records
organizedData$addBias(datasetNames = 'eBird')

}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-updateFormula"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-updateFormula}{}}}
\subsection{Method \code{updateFormula()}}{
Function used to update the formula for a selected observation model. The function is designed to work similarly to the generic \code{update} formula, and should be used to thin terms out of a process from the full model specified in \code{\link{intModel}}. The function also allows the user to add their own formula to the model, such that they can include non-linear components in the model. The function can also be used to print out the formula for a process by not specifying the \code{Formula} or \code{newFormula} arguments.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$updateFormula(
  datasetName = NULL,
  Formula,
  newFormula,
  processLevel = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetName}}{Name of the dataset (class \code{character}) for which the formula needs to be changed.}

\item{\code{Formula}}{An updated formula to give to the process. The syntax provided for the formula in this argument should be identical to the formula specification as in base \strong{R}. Should be used to thin terms out of a formula but could be used to add terms as well. If adding new terms not specified in \code{intModel}, remember to add the associated component using \code{.$changeComponents} as well.}

\item{\code{newFormula}}{Completely change the formula for a process -- primarily used to add non-linear components into the formula. Note: all terms need to be correctly specified here.}

\item{\code{processLevel}}{Logical argument: if \code{TRUE} changes the formulas for all of the processes in a dataset. Defaults to \code{FALSE}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An updated formula.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, 
                           responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Remove Forest from eBird
 organizedData$updateFormula(datasetName = 'eBird', Formula = ~ . - Forest)
 
 #Add some scaling to Forest for Parks
 organizedData$updateFormula('Parks', newFormula = ~ I(. +(Forest+1e-6)*scaling))
 
 #Now dd scaling to components
 organizedData$changeComponents(addComponent = 'scaling') 
 
}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-changeComponents"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-changeComponents}{}}}
\subsection{Method \code{changeComponents()}}{
Function to add and specify custom components to model, which are required by \pkg{inlabru}. The main purpose of the function is to re-specify or completely change components already in the model, however the user can also add completely new components to the model as well. In this case, the components need to be added to the correct formulas in the model using the \code{.$updateFormula} function. If \code{addComponent} and \code{removeComponent} are both missing, the function will print out the components to be supplied to \pkg{inlabru}'s \code{\link[inlabru]{bru}} function.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$changeComponents(addComponent, removeComponent, print = TRUE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{addComponent}}{Component to add to the integrated model. Note that if the user is re-specifying a component already present in the model, they do not need to remove the old component using \code{removeComponent}.}

\item{\code{removeComponent}}{Component (or just the name of a component) present in the model which should be removed.}

\item{\code{print}}{Logical: should the updated components be printed. Defaults to \code{TRUE}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An updated components list.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{

 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh, 
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Remove Forest from components
 organizedData$changeComponents(removeComponent = 'Forest')

}

}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-priorsFixed"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-priorsFixed}{}}}
\subsection{Method \code{priorsFixed()}}{
Function to change priors for the fixed (and possibly random) effects of the model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$priorsFixed(
  Effect,
  datasetName = NULL,
  mean.linear = 0,
  prec.linear = 0.001
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{Effect}}{Name of the fixed effect covariate to change the prior for. Can take on \code{'intercept'}, which will change the specification for an intercept (specified by one of \code{species} or \code{datasetName}).}

\item{\code{datasetName}}{Name of the dataset for which the prior of the intercept should change (if fixedEffect = 'intercept'). Defaults to \code{NULL} which will change the prior effect of the intercepts for all the datasets in the model.}

\item{\code{mean.linear}}{Mean value for the prior of the fixed effect. Defaults to \code{0}.}

\item{\code{prec.linear}}{Precision value for the prior of the fixed effect. Defaults to \code{0.001}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
New priors for the fixed effects.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Add prior to Forest
 organizedData$priorsFixed(Effect = 'Forest', mean.linear = 2, prec.linear = 0.1)

}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-specifySpatial"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-specifySpatial}{}}}
\subsection{Method \code{specifySpatial()}}{
Function to specify random fields in the model using penalizing complexity (PC) priors for the parameters.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$specifySpatial(
  sharedSpatial = FALSE,
  datasetName,
  Bias,
  PC = TRUE,
  Remove = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{sharedSpatial}}{Logical: specify the shared spatial field in the model. Requires \code{pointsSpatial == 'shared'} in \code{\link{intModel}}. Defaults to \code{FALSE}.}

\item{\code{datasetName}}{Name of which of the datasets' spatial fields to be specified. Requires \code{pointsSpatial = 'individual'} in \code{\link{intModel}}.}

\item{\code{Bias}}{Logical: specify the spatial field for the bias effect. If seperate fields are specified for different fields, the argument may be the name of the dataset for which the bias field to be specified.}

\item{\code{PC}}{Logical: should the Matern model be specified with pc priors. Defaults to \code{TRUE}, which uses \code{\link[INLA]{inla.spde2.pcmatern}} to specify the model; otherwise uses \code{\link[INLA]{inla.spde2.matern}}.}

\item{\code{Remove}}{Logical: should the chosen spatial field be removed. Requires one of \code{sharedSpatial}, \code{species}, \code{mark} or \code{bias} to be non-missing, which chooses which field to remove.}

\item{\code{...}}{Additional arguments used by \pkg{INLA}'s \code{\link[INLA]{inla.spde2.pcmatern}} or \code{\link[INLA]{inla.spde2.matern}} function, dependent on the value of \code{PC}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new model for the spatial effects.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Specify the shared spatial field
 organizedData$specifySpatial(sharedSpatial = TRUE,
                       prior.range = c(1,0.001),
                       prior.sigma = c(1,0.001))

} 
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-changeLink"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-changeLink}{}}}
\subsection{Method \code{changeLink()}}{
Function used to change the link function for a given process.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$changeLink(datasetName, Link)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetName}}{Name of the dataset for which the link function needs to be changed.}

\item{\code{Link}}{Name of the link function to add to the process. If missing, will print the link function of the specified dataset.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new link function for a process.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh, 
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present')

 #Specify the shared spatial field
 organizedData$changeLink('Parks', 'logit')
 
 
} 
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-spatialBlock"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-spatialBlock}{}}}
\subsection{Method \code{spatialBlock()}}{
Function to spatially block the datasets, which will then be used for model cross-validation with \code{\link{blockedCV}}. See the \code{\link[blockCV]{spatialBlock}} function from \pkg{blockCV} for how the spatial blocking works and for further details on the function's arguments.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$spatialBlock(k, rows_cols, plot = FALSE, seed = 1234, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{k}}{Integer value reflecting the number of folds to use.}

\item{\code{rows_cols}}{Integer value by which the area is divided into longitudinal and latitudinal bins.}

\item{\code{plot}}{Plot the cross-validation folds as well as the points across the boundary. Defaults to \code{FALSE}.}

\item{\code{seed}}{Seed used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}} to make the spatial blocking reproducible across different models. Defaults to \code{1234}.}

\item{\code{...}}{Additional arguments used by \pkg{blockCV}'s \code{\link[blockCV]{spatialBlock}}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
If \code{plot = TRUE}, a plot of the grid.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{ if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 Forest <- terra::rast(
 system.file(
 'extdata/SolitaryTinamouCovariates.tif', 
 package = "PointedSDMs"))$Forest
 
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           spatialCovariates = Forest,
                           Projection = proj, responsePA = 'Present',
                           pointsSpatial = 'individual')

 #Specify the spatial block
 organizedData$spatialBlock(k = 2, rows = 2, cols = 1, plot = FALSE)
 
} 

}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-addSamplers"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-addSamplers}{}}}
\subsection{Method \code{addSamplers()}}{
Function to add an integration domain for the PO datasets.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$addSamplers(datasetName, Samplers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetName}}{Name of the dataset for the samplers.}

\item{\code{Samplers}}{A \code{Spatial*} object representing the integration domain.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
New samplers for a process.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           Projection = proj, 
                           responsePA = 'Present')
 
#Add integration domain for the eBird records
organizedData$addSamplers(datasetName = 'eBird', Samplers = SolitaryTinamou$region)

}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-specifyRandom"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-specifyRandom}{}}}
\subsection{Method \code{specifyRandom()}}{
Function to specify the models and priors for the random effects included in the model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$specifyRandom(
  temporalModel = list(model = "ar1"),
  copyModel = list(beta = list(fixed = FALSE)),
  copyBias = list(beta = list(fixed = FALSE))
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{temporalModel}}{List of model specifications given to the control.group argument in the time effect component. Defaults to \code{list(model = 'ar1')}; see \code{\link[INLA]{control.group}} from the \pkg{INLA} package for more details.}

\item{\code{copyModel}}{List of model specifications given to the hyper parameters for the \code{"copy"} model. Defaults to \code{list(beta = list(fixed = FALSE))}.}

\item{\code{copyBias}}{List of model specifications given to the hyper parameters for the \code{"copy"} bias model. Defaults to \code{list(beta = list(fixed = FALSE))}.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An updated component list.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\dontrun{
 if (requireNamespace('INLA')) {
   
 #Get Data
 data("SolitaryTinamou")
 proj <- "+proj=longlat +ellps=WGS84"
 data <- SolitaryTinamou$datasets
 mesh <- SolitaryTinamou$mesh
 mesh$crs <- proj
 
 #Set model up
 organizedData <- startISDM(data, Mesh = mesh,
                           Projection = proj, 
                           responsePA = 'Present',
                           pointsSpatial = copy)
 
#Add integration domain for the eBird records
organizedData$specifyRandom(copyModel =  list(beta = list(fixed = TRUE)))

}
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-new"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$new(
  data,
  projection,
  Inlamesh,
  initialnames,
  responsecounts,
  responsepa,
  pointcovariates,
  trialspa,
  spatial,
  intercepts,
  spatialcovariates,
  boundary,
  ips,
  temporal,
  temporalmodel,
  offset,
  copymodel,
  formulas
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-specifyISDM-samplingBias"></a>}}
\if{latex}{\out{\hypertarget{method-specifyISDM-samplingBias}{}}}
\subsection{Method \code{samplingBias()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{specifyISDM$samplingBias(datasetName, Samplers)}\if{html}{\out{</div>}}
}

}
}
