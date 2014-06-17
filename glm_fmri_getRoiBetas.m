function B = glm_fmri_getRoiBetas(roiNii,betaDir,betaFStr,omitOutliers)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get voxelwise beta estimates within an roi mask and average over them


% INPUTS:
%     roiNii  - either nii or filepath to nii roi file. the file shoudl
%     contain a zero or nan at each voxel not included in the roi and a
%     non-zero value for all voxels within the roi. In most cases this will
%     be a binary mask of zeros and ones. The non-zero values will be used
%     for a weighted-mean across the roi voxels.
%
%     betaDir - string specifying the directory containing the nifti beta
%               files
%     betaFStr- file identifier string specifying the files to get betas
%               from, e.g., 'shock','sn_corr',etc. Will also use subject
%               strings to id files to use.
%    omitOutliers (optional) - if set to 1, voxel betas with an
%               abs(zscore)>3 will be omitted from the roi beta estimation.
%               Default is 0 (to not omit outliers)


% OUTPUTS:
%     B - n x p matrix with n=# of subjects and p=# of betas


% kjh, 09-Dec-2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get roi mask


% roi nifti file
if ischar(roiNii)
    roiNii = readFileNifti(roiNii);
end


% dont omit outliers unless omitOutliers = 1
if ~exist('omitOutliers') || isempty(omitOutliers)
    omitOutliers = 0;
end


% get all list of all beta files to read
a = dir([betaDir,'/*',betaFStr,'*']);


% subjects loop
for s = 1:length(a)
    
    nii = readFileNifti(fullfile(betaDir,a(s).name));
    
    % if betas are from spm_hrf files, use only the first vol
    if strfind(betaDir,'spm_hrf')
        nii.data = nii.data(:,:,:,1);
    end
    
    % extract betas from roi voxels
    vox_betas=reshape(nii.data,prod(nii.dim(1:3)),[]);
    vox_betas=vox_betas(roiNii.data~=0,:);
    
    
    % for a weighted mean; if roi mask is binary, then its a simple mean
    w = roiNii.data(roiNii.data~=0);
    w = repmat(w,1,size(vox_betas,2));
    
    
    if omitOutliers
        
        oidx=find(abs(zscore(vox_betas,0,2))>3);
        vox_betas(oidx) = nan;
        w(oidx) = nan;
        
    end
    
    % average across non-zero (roi) voxels. The mean is weighted by the
    % values in roiNii.data
    B(s,1:size(vox_betas,2))= (nansum(w.*vox_betas)./nansum(w));
    
    
end


B(B==0)=nan;  % for per trial betas


