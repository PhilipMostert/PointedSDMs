# PointedSDMs 2.1.2

#### Changes and fixes since previous version:

-   Various fixes to `predict` for the multi-species model.
-   Documentation updates.
-   Fixes to the two cross-validation functions. Most notably in `blockedCV`. A new method has been added which allows the user to compare different combinations of datasets in the model by leaving datasets out of a train model, predicting on a dataset, and using the prediction as an offset in a testing model.\
