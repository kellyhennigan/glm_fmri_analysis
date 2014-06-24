% script to calculate seed ROI-voxel beta series correlations

% loops through subjects to calculate those correlations then conducts
% tests for differences in correlations between conditions (defined by
% stims).

% Tests are conducted on Fisher's r-to-z transformed to correlation
% coefficients.

% Will output a nifti file with the test statistics.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

clear all
close all

inDir = ['/Users/Kelly/ShockAwe/data/results_per_trial'];  % relative to subject directory
outDir = ['/Users/Kelly/ShockAwe/data/ttests_per_trial'];  % relative to subject directory

saveStatMap = 1; % if 1/true, save out stat maps as nifti files

maskFile = '/Users/Kelly/ShockAwe/data/ROIs_tlrc/group_mask.nii.gz';

stims = {'shock','neutral'};

N = 18;

seedStr = 'peakLHab'; % string to use for out files

% provide string specifying the directory where seed beta series are
% saved
seedBDir = '/Users/Kelly/ShockAwe/data/betas/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute beta series correlations


mask = readFileNifti(maskFile); % get mask file
mask.data = double(mask.data);
idx = find(mask.data);


for c = 1:length(stims)
    
    % load seed betas
    seedBs = dlmread([seedBDir seedStr '_' stims{c} '_per_trial_betas'])';
    
    % load voxel betas
    cd(inDir);
    f = dir(['sa*',stims{c},'*nii*']);
    N = length(f);
    
    for i = 1:N
          
        % get voxel beta series for this subject
        nii = readFileNifti(f(i).name);
        subData = reshape(nii.data,prod(nii.dim(1:3)),[]); % each volume of voxels is in a row now
        subData = subData(idx,:)'; % get only voxels within the mask
        
        B2 = permute(subData,[1 3 2]); % reshape so voxel beta series are in the 3rd dim for calcBSCorr()
        
        % get this subject's seed ROI beta series
        B1 = seedBs(1:size(B2,1),i);
      
        
        [r(i,c,:),Z(i,c,:)] =calcBSCorr(B1,B2);
        
        clear B1 B2
        
    end % subjects
    
    
end % stims


%% tests

switch numel(stims)
    
    case 1      % H0: B1 = 0
        
        [~,p,~,stats] = ttest(Z);
        testStat = stats.tstat;
        df = mode(stats.df); % take the mode in case of random NaN voxels
        
        % mean Z across subjects for every voxel
        meanStatMap = mask.data;
        meanStatMap(idx) = mean(Z,1);
        
        descrip = ['one sample t-test; df(' num2str(df) '); vol1=t-stat; vol2=p'];
        
    case 2  % H0: B1 = B2
        
        [~,p,~,stats] = ttest(Z(:,1,:),Z(:,2,:));
        testStat = stats.tstat;
        df = mode(stats.df); % take the mode in case of random NaN voxels
        
        % mean Z across subjects for every voxel
        meanStatMap = mask.data;
        meanStatMap(idx) = nanmean(Z(:,1,:)-Z(:,2,:),1);
        
        descrip = ['paired sample t-test; df(' num2str(df) '); vol1=t-stat; vol2=p'];
        
    case 3 % H0: B1 = B2 = B3
        
        % repeated measures ANOVA
        % how can I do this without using a for loop?? :/
        for v = 1:size(Z,3)
            vox_betas = Z(:,:,v);
            [p(v),table]=anova_rm(Z(:,:,v),'off');
            testStat(v) = table{2,5};  % F-statistic
            df1(v)=table{2,3}; df2(v)=table{4,3}; % degrees of freedom (between conds, error)
        end
        df(1) = mode(df1); df(2) = mode(df2);
        
        meanStatMap = []; % not sure what makes sense for a meanStatMap here
        
        descrip = ['rep_meas_ANOVA; df(' num2str(df(1)) ',' num2str(df(2)) '); vol1=Fstat; vol2=p'];
end

% test-stat (F or t)
testStatMap = mask.data;
testStatMap(idx) = testStat;

% p-value map
pMap = mask.data;
pMap(idx) = p;


if saveStatMap
    cd(outDir);
    outNii = makeGlmNifti(mask,[seedStr '_' stims{1},'-',stims{2}],descrip,testStatMap,pMap);
    writeFileNifti(outNii);
end



