function roiNiftiToMat(niiRoiFile, outName)
%
% converts a nifti roi file to a .mat file 
% of the same format as those generated 
% with dtiFiberUI and saves it in the same directory
% 
% uses the dtiNewRoi function to make the Roi
%
% inputs: niiRoiFile - file name of nii file (e.g., 'ROI.nii.gz')
%         
% outputs: roi .mat file of the same name in the same directory
% 
% note: if outName argument is omitted, then the new .mat roi file will
% take the same name as the nifti version
%
% kjh 4/2011
%
% I took code from dtiConvertFreeSurferRoiToMat
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

[path, fname] = fileparts(niiRoiFile); % if only filename is given than path returns ''

strIndx = strfind(niiRoiFile,'.nii');
if strIndx
    fname = fname(1:strIndx-1);
end

    
if (~exist('outName','var') || isempty(outName))
  outName = fname;
end

% load nifti ROI
niiRoi = readFileNifti(niiRoiFile);

% if the roi file doesn't load, throw an error
if isempty(niiRoi.data)
    error('roi file isn''t loading');
end

% get ROI coords
indx = find(niiRoi.data);
[i j k] = ind2sub(size(niiRoi.data), indx);

%convert to acpc coords
acpcCoords=round(mrAnatXformCoords(niiRoi.qto_xyz, [i j k]));


% create a new roi with mrDiffusion structure
roi=dtiNewRoi(outName,[],acpcCoords); % name, color, coords

dtiWriteRoi(roi,fullfile(path,[outName,'.mat']), [],'acpc');  % roi, filename, versionNum, coordinateSpace, xform
fprintf('\nwriting file %s\n',[outName,'.mat']);

