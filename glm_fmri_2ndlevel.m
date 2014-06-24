% % glm_fmri_2ndlevel_script
% --------------------------------
% usage: script to define group-level contrasts of beta estimates and test
% them

% NOTES: to do: add functionality for testing linear contrasts, etc.

% kjh, June 2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

irf = 'spm_hrf';

mainDir = '/Users/Kelly/ShockAwe/data/';

inDir = [mainDir 'results_' irf];  % relative to subject directory
outDir = [mainDir 'ttests_' irf];  % relative to subject directory
saveStatMap = 1; % if 1/true, save out stat maps as nifti files

maskFile = [mainDir 'ROIs_tlrc/group_mask.nii.gz'];


% these strings will be used to identify nifti files to load from inDir
stims = {'juice','neutral','shock'};
volIdx = 1;   % volume index specifying which volume contains stat to test
% this should be either an integer or the string 'peak' will find the
% volume with the highest beta estimate


%% load nifti files

% get mask
mask = readFileNifti(maskFile);
mask.data = double(mask.data);
idx = find(mask.data);


% get stat files
cd(inDir);

for j = 1:numel(stims)
    
    f = dir(['sa*',stims{j},'*.nii*']);
    N = length(f);
    fprintf(['\n\nloading ' stims{j} ' files...']);
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
    fprintf('done.');
    
    meanStatMap{j} = mask.data;
    meanStatMap{j}(idx) = mean(d{j});
    
end  % stims



%% tests


switch numel(stims)
    
    case 1  % H0: B1 = 0
        
        [~,p,~,stats] = ttest(d{1});
        testStat = stats.tstat;
        df = stats.df;
        
        meanStatMap = mask.data;
        meanStatMap(idx) = mean(d{1});
        
        
    case 2  % H0: B1 = B2
        
        [~,p,~,stats] = ttest(d{1},d{2}); % paired ttest
        testStat = stats.tstat;
        df = stats.df;
        
        meanStat = mean(d{1}-d{2});
        meanStatMap = mask.data;
        meanStatMap(idx) = meanStat;
        
        
    case 3 % H0: B1 = B2 = B3
        
        % there's gotta be a way to do this without a voxelwise for loop!?
        for j = 1:size(d{1},2)
            % repeated measures ANOVA
            [this_p,table]=anova_rm([d{1}(:,j),d{2}(:,j),d{3}(:,j)],'off');
            p(j) = this_p(1);
            testStat(j) = table{2,5};  % F-statistic
        end
        df(1)=table{2,3}; % degrees of freedom (between conds, error)
        df(2)=table{4,3};
        
        meanStatMap = []; % not sure what makes sense for a meanStatMap here
        
end

% make test stat (t or F) maps and p maps
testStatMap =  mask.data;
testStatMap(idx) = testStat;

pMap = mask.data;
pMap(idx) = p;



%% save out test results

% save out nifti file & convert to afni file
if(saveStatMap)
    
    cd(outDir);
    
    if numel(stims) <= 2  % t-test
        outNameStr = cat(2,stims{:});
        descrip = ['one sample or paired ttest;  df(' num2str(df) ')'];
        outNii = makeGlmNifti(mask,outNameStr,descrip,meanStatMap,testStatMap);
    else
        outNameStr = [outNameStr '_F_p']; % repeated measures one-way ANOVA
        descrip = ['rep_meas_ANOVA; df(' num2str(df(1)) ',' num2str(df(2)) ')'];
        outNii = makeGlmNifti(mask,outNameStr,descrip,testStatMap,pMap);
    end
    
    writeFileNifti(outNii);
    
end

