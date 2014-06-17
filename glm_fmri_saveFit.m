% fit the model to pre-processed fmri data and save out the results


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

subjects = getSASubjects('fmri');

irf='per_trial';
% irf = 'spm_hrf';
irf_param_str = '_6_16_1_1_10_0_32';
% irf = 'tent';
% irf_param_str = '_0_10_6';

gSpace = 'tlrc';  % '' for native space

matName = ['design_mats/glm_',irf,irf_param_str,'wTRs.mat'];

funcFile = ['afni/all_scaled_s5_',gSpace,'.nii.gz'];

maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/group_mask.nii.gz';

outDir = ['/home/kelly/ShockAwe/data/results_',irf,'_s5'];  % relative to subject directory

stims = {'juice','neutral','shock'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subject loop

mask=readFileNifti(maskFile);

for s =2:18
    
    subject = subjects{s};
    
    expPaths = getSAPaths(subject);
    cd(expPaths.subj);
    
    % get design matrix & data
    load(matName);
    func = readFileNifti(funcFile);
    
    % remove TRs to censor from the data & the associated rows from the design matrix
    [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx);
    
    % fit model to data
    stats = glm_fmri_fit_vol(data,X,regIndx,mask.data);
    
    
    cd(outDir);
    
    for c = 1:length(stims)
        outName = [subject,'_',stims{c},'_betas'];
        outNii = makeGlmNifti(mask,outName,stats.B(:,:,:,[strmatch(stims{c},regLabels)]));
        writeFileNifti(outNii); 
    end
    
    
    clear expPaths func X
    
    fprintf(['\n\nfinished subject ', subject,'\n']);
    
    
end % subjects
