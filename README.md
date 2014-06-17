These scripts and sample data go through the steps for conducting a beta series correlation analysis, a measure of functional connectivity. 

Steps: 

1) generate regressor time series
2) construct glm
3) fit model to sample data
4) correlate ROI-voxel beta series to create a voxelwise map of correlation values

once this is done for each subject, transform these r values (Fisher's r-to-z transform) and then you can do voxelwise t-tests to test whether the correlations during condition A and B differ. 