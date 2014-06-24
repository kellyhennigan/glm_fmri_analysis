function outNii = makeGlmNifti(refNii, fname, descrip, varargin)
% --------------------------------
% usage: creates a new nifti struct that saves voxelwise stat maps, etc.


% INPUT:
%   refNii - reference nifti to get all header info (xform, voxel
%            dimensions, etc.)
%   fname - string specifying the file name for the new out nifti
%           descrip (optional) - string giving additional info about this
%           data, e.g., 'vol1 is Fstat; vol2 is p-stat; df(2,17)
%   varargin - 3d data matrices to include in out nii; must be the same
%           dimensions as the refNii

% OUTPUT:
%   outNii - nifti struct containing info/data specified by input

% NOTES:
%
% kelly 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outNii = refNii; % get all header info from refNii

% if fname is not defined, call it outNii
if notDefined('fname')
    fname = 'outNii.nii.gz';
end

% if fname doesn't have nii file extension, had it
suffix = strfind(fname,'.nii');
if isempty(suffix)
    fname = [fname,'.nii.gz'];
end
outNii.fname = fname;


% fill in data w/input volumes
outNii.data = [];
if ~isempty(varargin)
    for h = 1:length(varargin)
        outNii.data = cat(4,outNii.data, varargin{h});
    end
end

% add description of the data if provided
outNii.descrip = descrip;

outNii.dim = size(outNii.data);

