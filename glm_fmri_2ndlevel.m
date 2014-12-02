% % glm_fmri_2ndlevel_script
% --------------------------------
% usage: script to define group-level contrasts of beta estimates and test
% them

% NOTES: to do: add functionality for testing linear contrasts, etc.

% kjh, June 2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

clear all
close all

irf = 'tent';

mainDir = '/Users/Kelly/ShockAwe/data/';

inDir = [mainDir 'results_' irf '_new'];  % directory containing subject beta maps contains bet
outDir = [mainDir 'ttests_' irf '_new'];  % directory to save out results to

saveStatMap = 1; % if 1/true, save out stat maps as nifti files

maskFile = [mainDir 'ROIs_tlrc/group_mask.nii.gz'];


% these strings will be used to identify nifti files to load from inDir
stims = {'shock','neutral'};
volIdx = 2;   % volume index specifying which volume contains stat to test
% this should be either an integer or the string 'peak' will find the
% volume with the highest beta estimate


%% load nifti files

% get mask
mask = readFileNifti(maskFile);
mask.data = double(mask.data);
idx = find(mask.data);


% get stat files
cd(inDir);

fprintf('\n\nloading stim files...');

for j = 1:numel(stims)
    
    
    f = dir(['sa*',stims{j},'*.nii*']);
    N = length(f);
    
    for i = 1:N
        nii = readFileNifti(f(i).name);
        subData = reshape(nii.data,prod(nii.dim(1:3)),[]); % each volume of voxels is in a row now
        subData = subData(idx,:);  % get only voxels within the mask
        
        if(strcmp(volIdx,'peak'))
            [peakB, peakBIdx]=max(subData(:,3:5),[],2);  % only look for peak between 4-8 sec after stim onset
            all_peakBIdx(:,i) = peakBIdx+2; % add 2 to index to make up for not using the first 2 regressors at time 0 or 2
            d{j}(i,:) = peakB;
        else
            d{j}(i,:) = mean(subData(:,volIdx),2); % take mean over columns in case there's more than 1 volIdx
        end
        
    end  % subjects
    
    meanStatMap{j} = mask.data;
    meanStatMap{j}(idx) = mean(d{j});
    
end  % stims

fprintf('done.');


%% tests

fprintf('\n\nconducting tests...');

switch numel(stims)
    
    case 1  % H0: B1 = 0
        
        [~,p,~,stats] = ttest(d{1});
        testStat = stats.tstat;
        df = mode(stats.df);
        
        outNameStr = [stims{1} '_T_p'];
        descrip = ['one sample t-test;  df(' num2str(df) ')'];
        
        
    case 2  % H0: B1 = B2
        
        [~,p,~,stats] = ttest(d{1},d{2}); % paired ttest
        testStat = stats.tstat;
        df = mode(stats.df);
        
        outNameStr = [stims{1} '-' stims{2} '_T_p'];
        descrip = ['paired t-test;  df(' num2str(df) ')'];
        
        
    case 3 % H0: B1 = B2 = B3
        
        % there's gotta be a way to do this without a voxelwise for loop!?
        for j = 1:size(d{1},2)
            % repeated measures ANOVA
            [this_p,table]=anova_rm([d{1}(:,j),d{2}(:,j),d{3}(:,j)],'off');
            p(j) = this_p(1);
            testStat(j) = table{2,5};  % F-statistic
            df1(j)=table{2,3}; df2(j)=table{4,3}; % degrees of freedom (between conds, error)
        end
        df(1) = mode(df1); df(2) = mode(df2);
        
        outNameStr = [stims{1} '-' stims{2} '-' stims{3} '_F_p'];
        descrip = ['rep_meas_ANOVA; df(' num2str(df(1)) ',' num2str(df(2)) '); vol1=Fstat; vol2=p'];
        
end


% make test stat (t or F) maps and p maps
testStatMap =  mask.data;
testStatMap(idx) = testStat;

pMap = mask.data;
pMap(idx) = p;


fprintf('done.');


%% save out test results

% save out nifti file & convert to afni file
if(saveStatMap)
    
    fprintf('\n\nsaving results...');
    
    cd(outDir);
    
    outNii = makeGlmNifti(mask,outNameStr,descrip,testStatMap,pMap);
    
    writeFileNifti(outNii);
    
    fprintf('done.');
    
end

