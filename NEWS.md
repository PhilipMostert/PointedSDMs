# PointedSDMs 1.2.0

#### Fixes since previous version:

-   New vignette (*Marked_Point_Process)* which shows off the marked point process component of the model. This vignette also has a new dataset (*Koala*) which includes observations of blue gum across Phillip island (Australia) with various marks describing each sighting.
-   Fixed issue regarding adding *pointCovariates* to *PO* data in `intModel`.
-   Added `marks` to `predict` for objects of class `bruSDM`, which allows the user to make predictions for the marks in the model.
-   Added temporal model to the marks in `intModel`.
-   Added *Copy* in `.$specifySpatial`, which allows the user to use *INLA*'s copy feature for the marks and bias spatial random fields.
-   Various documentation changes.
