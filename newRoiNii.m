function roiNii = newRoiNii(templateNii, coordIdx, outName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to make a new nifti roi mask file


% INPUTS:
%     templateNii - a nifti file to use as a template that has the same
%                   dimensions, xforms, etc. that the new roi should have 
%     coordIdx -  image coordinates of the voxels to include in the ROI. Can
%     be an index or [i,j,k] subscript.
%     outName - string that gives the file name for the new ROI nifti file.
%     If left blank, the returned nifti file will be called 'ROI.nii.gz'.


% OUTPUTS:
%     roiNii - new nifti ROI file 


% kjh, Aug 2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % define new roi nifti from template nifti 
roiNii = templateNii; 
roiNii.data = zeros(size(roiNii.data(:,:,:,1)));
roiNii.descrip = 'roi mask';


% set up coordIdx
if notDefined('coordIdx')
    coordIdx = [];
end

% get roi coordinate voxel indices
if any(size(coordIdx) == 3)
    if size(coordIdx,2)~=3
        coordIdx = coordIdx';
    end
    coordIdx = sub2ind(size(roiNii.data),coordIdx(:,1),coordIdx(:,2),coordIdx(:,3));
end
roiNii.data(coordIdx) = 1;


% set up outName
if notDefined('outName')
    outName = 'ROI.nii.gz';
end

if isempty(strfind(outName,'.nii'))
    outName = [outName '.nii.gz'];
end

roiNii.fname = outName;


    