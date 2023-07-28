# PointedSDMs 1.3.0

#### Changes and fixes since previous version:

-   Migration away from the R packages *sp* and *Raster* towards the more up-to-date *sf* and *terra*. The default class for species occurrence data should now be *sf* objects, and the default for spatial covariates should be a *spatRaster* object. This shift has not changed the fundamentals of the package, nor has it changed any of the function's arguments.
-   Changed all of the vignettes and their underlying data as a result of this shift.
-   Various internal code changes as a result of changes in *inlabru* version *2.8.0*.
-   Fixed the predict function for when `predictor = TRUE.`
-   Fixed the `.spatialBlock()` function to work with *blockCV*'s `cv_spatial`.
-   Various spelling and grammar fixes in documentation.
