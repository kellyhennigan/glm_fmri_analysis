function outNii = makeGlmNifti(refNii, fname, varargin)

% creates a new nifti struct named fname w/ refNii header info and 
% data given in the last arguments concatenated along the fourth dimension 
% in the order they are given.

% kelly 2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outNii = refNii;

suffix = strfind('.nii', fname);
if isempty(suffix)
    fname = [fname,'.nii.gz'];
end

outNii.fname = fname;
outNii.data = [];

% fill in data
    if ~isempty(varargin)
        for h = 1:length(varargin)
           
            outNii.data = cat(4,outNii.data, varargin{h});
        
        end
        
    end
    
    outNii.dim = size(outNii.data);
    
    