% make a glm design matrix and save it out as a mat file along with
% regressor labels and an index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files, params, etc

subj = '18';

[~,cb]=getSA2Subjects(subj);

context_runs = {'base','base','base','stress','stress','stress'}; % baseline or stress context

runs = 1:6; % scan runs to include

pa=getSA2Paths(subj);

nVols = 326; % vector w/corresponding # of vols/run, or one scalar to use the same value for all runs

irf = 'spm_hrf'; 
% irfParamsStr = '';

outMatStr = ['glm_',irf, irfParamsStr];
% outMatName = 'glm_habPPI.mat';

%%%%% Regressors of interest %%%%%

stims = {'gain+1','gain+PE','gain0','gain-PE','loss-1','loss-PE','loss0','loss+PE',...
    'contextevent','shock','cuepair1','cuepair2'};


%%%%% Baseline/Nuisance regressors %%%%%

motionRegsFiles = {'vr_run.1D','vr_run2.1D','vr_run3.1D','vr_run4.1D','vr_run5.1D','vr_run6.1D'};

% number of polynomial baseline regressors to include for each scan run
nPolyRegs = 5; % should be at least 1, +1 more for every 2.5 min of continuous scanning

censor_trs = []; % vols to censor bc of bad movement or otherwise



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do it
%  
    
    pa = getSA2Paths(subj);
    
    %% Design matrix
    
    [X, regLabels, regIndx] = glm_fmri_mat(pa, subject, runs, stims, ...
        motionRegsFiles, nPolyRegs, censor_first_trs);
    
    %% save it
    
    cd(expPaths.design_mats); 
    save([outMatName,'.mat'],'X','regLabels','regIndx');
    
    fprintf(['\n\nfinished subject', subject]);
    
end % subjects