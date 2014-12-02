% script for fitting a glm to subjects' pre-processed fmri data and
% saving the results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

subjects = getSASubjects('fmri');

% irf='per_trial';
% irf = 'spm_hrf';
% irf_param_str = '_6_16_1_1_6_0_32';
irf = 'tent';
irf_param_str = '_0_10_6';

gSpace = 'tlrc';  % '' for native space

matName = ['glm_',irf,irf_param_str '.mat'];

funcFile = ['all_scaled_s' gSpace '.nii.gz'];

maskFile = '/Users/Kelly/ShockAwe/data/ROIs_tlrc/new/llatSN_F.nii.gz';

outDir = ['/Users/Kelly/ShockAwe/data/results_' irf '_noTRs'];  % relative to subject directory

stims = {'juice','neutral','shock'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subject loop

mask=readFileNifti(maskFile);

for s =1:18

    
    subject = subjects{s};
    
    expPaths = getSAPaths(subject);
    cd(expPaths.subj);
    
    % get design matrix
    matPath = fullfile(expPaths.design_mats,matName);
    load(matPath);
   
    
    % get data
    func = readFileNifti(funcFile);
    
    % remove TRs to censor from the data & the associated rows from the design matrix
    [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx);
    
    % fit model to data
    stats = glm_fmri_fit_vol(data,X,regIndx,mask.data);
    
    cd(outDir);
    
    for c = 1:length(stims)
        

        out_descrip = ['glm file: ' matName];
       
        outName = [subject '_' stims{c} '_betas'];
        outNii = makeGlmNifti(mask,outName,out_descrip,stats.B(:,:,:,[strmatch(stims{c},regLabels)]));
        writeFileNifti(outNii);
        
        
    end
    
    clear expPaths func X
    
    fprintf(['\n\nfinished subject ', subject,'\n']);
    
    
end % subjects
