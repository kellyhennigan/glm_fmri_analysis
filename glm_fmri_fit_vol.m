function stats = glm_fmri_fit_vol(data,X,regIdx,mask)
%
% this function fits a general linear model to a 3d volume of voxels w/
% time in the 4th dimension using least squares & returns the estimated
% parameters. It skips voxels where the mask (if given) is set to 0 or
% false.

% this function basically reshapes the input data so that each voxel is a
% column with time series as rows. It then calls glm_fmri_fit to estimate
% the model & reshapes data back to match the input dimensions.

%%%%%%%%%%%%%%%%%%%%%% INPUTS:

% data - 3d volume of voxel time series with time in the 4th dimension. 
% X - design matrix of regressors
% regIdx (optional) - column vector w/length = # of regressors with zeros 
%       for baseline regressors and 1,2, etc. for regressors of interest; 
% mask (optional) - 3D dataset of 0s and 1s matching data dimensions

%%%%%%%%%%%%%%%%%%%%%% OUTPUTS:

% a stats struct with the following fields:

% 'B'       - regressor coefficients estimated w/OLS
% 'seB'     - standard error of each beta coefficient
% 'tB'      - this returns the t or F-stat grouped by regIdx as appropriate
% 'pB'      - corresponding p-values for each t-stat
% 'Rsq'     - Rsquared for the full model not including the baseline
% 'err_ts'   - residual time series

% and if regIdx is is provided, also:

% 'Fstat'   - F-stat testing whether the full model is better than baseline
% 'pF'      - p value on F-statistic


% kelly May 2012

% Oct 2013, revised to perform all operations without loops

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check inputs & define some useful variables

dim = size(data);       % data dimensions
[nt,np] = size(X);      % # of modeled time points and parameters
nvox = prod(dim(1:3));  % total # of voxels within the 3d volume

% make sure size(X,1) equals size(data,4)
if nt ~= dim(4)
    error('the # of design matrix rows must equal the # of data timepoints (in the 4th dim).\n\n');
end


% if regIdx isnt given, set it to empty
if ~exist('regIdx','var')
    regIdx = [];
end


% if a mask is given, get an index of voxels to include
if exist('mask','var')
    if dim(1:3) ~= size(mask)
        error('data and mask dimensions must agree\n\n');
    end
    idx = find(mask~=0);
else  % if no mask is given, include all voxels
    idx = 1:nvox;
end


%% get to it

data = reshape(data,nvox,nt)'; % reshape data with time in rows & voxels as columns

stats = glm_fmri_fit(data(:,idx),X,regIdx); % fit model to masked data


%% reshape out stats to match the volume dimensions of the input data

fnames = fieldnames(stats);

for i = 1:numel(fnames)
    thisStat = zeros(size(getfield(stats,fnames{i}),1),nvox); % array of zeros
    thisStat(:,idx) = getfield(stats,fnames{i}); % fill in stats for modeled voxels 
    stats=setfield(stats,fnames{i},reshape(thisStat',[dim(1:3),size(thisStat,1)])); % put reshaped data in out struct
end
