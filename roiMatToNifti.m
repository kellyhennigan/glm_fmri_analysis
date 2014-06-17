function roiMatToNifti(matRoiFile, niiRefFile,outFName)
%
% converts a .mat roi to a nifti file 
% and saves it in the same directory
%
% note: assumes this is called from a subject's main directory and that the
% nii ref file is in the current dir
%
%
% inputs: matRoiFile - .mat roi file to convert 
%         refFile - reference file for header info
%         outFName (optional) - name of new nii file
%
% outputs: Roi nifti file in the current directory
%
% kjh 4/2011
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
    
ref = readFileNifti(niiRefFile);

load(matRoiFile); % loads a struct called 'roi'

imgCoords = round(mrAnatXformCoords(ref.qto_ijk, roi.coords));
coordIndx = sub2ind(ref.dim(1:3),imgCoords(:,1),imgCoords(:,2),imgCoords(:,3));

niiRoi = ref;
niiRoi.data = zeros(ref.dim(1:3));
niiRoi.data(coordIndx) = 1;

% new nii file name
if ~exist('outFName','var') || isempty(outFName)
    outFName = roi.name;
end
strIdx = strfind(outFName,'nii');
if isempty(strIdx)
    niiRoi.fname = [outFName,'.nii.gz'];
else
    niiRoi.fname = outFName;
end

% save new nifti ROI file
fprintf(['\n\n writing out ',niiRoi.fname,'\n\n']);
writeFileNifti(niiRoi);

