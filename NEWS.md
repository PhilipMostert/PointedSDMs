# PointedSDMs 1.3.2

#### Changes and fixes since previous version:

-   Added the argument *speciesIndependent* to `intModel`. This argument is logical and indicates if species effects should be made independently for species across datasets or not. Defaults to `FALSE`.
-   Added the argument *speciesEffects* to `intModel`. This argument takes a list with two named items: *randomIntercepts* and *Environmental*, indicating if species should have their own environmental and (random )intercept effects.
-   Added the option *shared* for *speciesSpatial* in `intModel`. This creates one random effect shared for all the species considered in the model.
-   Added *shareModel* to `$addBias`: allows the user to share bias fields across different datasets.
-   Changed the setup of the covariates: now completed in the dataset up stage rather than the modelling stage. Covariate data is now also added to the prediction data directly.
-   Updates to the *Solitary_tinamou* vignette to add new *pc* priors to the different models.
-   Removed *Raster* support for the spatial covariates
