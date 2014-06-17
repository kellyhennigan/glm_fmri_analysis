function B = glm_fmri_getRoiBetas2(roiFiles,roiStrs,irf)

% get the mean time series for a region of interest then fit a model to
% that time series


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files, params, etc

if ~iscell(roiFiles)
    roiFiles = {roiFiles};
end

if ~iscell(roiStrs)
    roiStrs = {roiStrs};
end


subjects = getSASubjects('fmri');

stims = {'juice','neutral','shock'}; % should correspond to regFiles

switch irf
    case {'spm_hrf','per_trial'}
        irf_param_str = '_6_16_1_1_10_0_32';
        
    case 'tent'
        irf_param_str = '_0_10_6';
end

design_mat = ['design_mats/glm_',irf,irf_param_str,'wTRs.mat'];

gSpace = 'tlrc';  % leave blank for native space, otherwise 'tlrc'


%% get roi mask if in group space

for r = 1:length(roiFiles)
    
    if (gSpace)
        fprintf(['\n\nworking on betas for roi ', roiFiles{r},'\n']);
        mask = readFileNifti(['/home/kelly/ShockAwe/data/ROIs_tlrc/',roiFiles{r}]);
        mask_idx = find(mask.data);
    end
    
    
    %% subject loop
    
    for s = 1:length(subjects)
        
        subject = subjects{s};
        
        fprintf(['\nWorking on subject ',subject,'...\n\n']);
        
        expPaths = getSAPaths(subject);
        cd(expPaths.subj);
        
        % get roi mask if in native space
        if ~(gSpace)
            mask=readFileNifti(fullfile(expPaths.ROIs,roiFile));
            mask_idx = find(mask.data);
        end
        
        if (gSpace)
            cd(['/home/kelly/ShockAwe/data/results_',irf]);
            
            
            switch irf
                
                %%%%%%%%%%%%%%%% spm_hrf
                
                case 'spm_hrf'
                    
                    if (gSpace)
                        cd(['/home/kelly/ShockAwe/data/results_',irf]);
                        for c = 1:length(stims)
                            nii = readFileNifti([subject,'_',stims{c},'_betaT.nii.gz']);
                            betaMap = nii.data(:,:,:,1);
                            B{c}(s,1) = mean(betaMap(mask_idx));
                        end
                    else
                        cd(expPaths.subj);
                        load(design_mat); % load design matrix
                        func = readFileNifti('afni/all_scaled_s.nii.gz'); % get func data
                        [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx); % censor first TRs
                        mean_ts = roi_mean_ts(data,mask.data); % get mean of ROI time series
                        stats = glm_fmri_fit(mean_ts,X,regIndx);
                        
                        for c = 1:length(stims)
                            idx = strmatch(stims{c},regLabels);
                            betas = stats.B(idx);
                            B{c}(s,1:length(betas))=betas;
                        end
                    end
                    
                    %%%%%%%%%%%%%%%% tent
                    
                case 'tent'
                    
                    if(gSpace)
                        cd(['/home/kelly/ShockAwe/data/results_',irf]);
                        for c = 1:length(stims)
                            nii = readFileNifti([subject,'_',stims{c},'_irf.nii.gz']);
                            betaMaps=reshape(nii.data,prod(nii.dim(1:3)),[]);
                            %                         betaMaps=betaMaps(:,1:2:size(betaMaps,2));
                            B{c}(s,:) = mean(betaMaps(mask_idx,:));
                        end
                    else
                        cd(expPaths.subj);
                        load(design_mat); % load design matrix
                        func = readFileNifti('afni/all_scaled_s.nii.gz'); % get func data
                        [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx); % censor first TRs
                        mean_ts = roi_mean_ts(data,mask.data); % get mean of ROI time series
                        stats = glm_fmri_stats(mean_ts,X,regIndx);
                        
                        for c = 1:length(stims)
                            idx = strmatch(stims{c},regLabels);
                            betas = stats.B(idx);
                            B{c}(s,1:length(betas))=betas;
                        end
                        
                    end
                    
                    %%%%%%%%%%%%%%%% per_trial
                    
                case 'per_trial'
                    
                    cd(expPaths.subj);
                    load(design_mat); % load design matrix
                    func = readFileNifti(['afni/all_scaled_s',gSpace,'.nii.gz']); % get func data
                    [data,X,regLabels,regIndx] = glm_fmri_censorTRs(func.data,X,regLabels,regIndx); % censor first TRs
                    mean_ts = roi_mean_ts(data,mask.data); % get mean of ROI time series
                    stats = glm_fmri_stats(mean_ts,X,regIndx);
                    
                    for c = 1:length(stims)
                        idx = strmatch(stims{c},regLabels);
                        betas = stats.B(idx);
                        B{c}(s,1:22)=nan;
                        B{c}(s,1:length(betas))=betas;
                    end
                    
            end % irf
            
        end % subjects
        
        
        B_all = [B{:}];
        alpha = .05;
        if (plotFigs)
            p_vals=testBetas(B)
            h(f)=glm_fmri_plotBetas(B,irf);
            title(gca,roiStr)
            %     h=glm_fmri_plotBetas(B,irf,p_vals);
        end
        
        
        cd(outDir);
        for c = 1:length(stims)
            dlmwrite(outFNames{c},B{c});
        end
    end
    
end

% B2=B;
