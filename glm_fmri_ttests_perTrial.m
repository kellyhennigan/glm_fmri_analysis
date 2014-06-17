% % glm_fmri_ttests for beta series correlations

% give a seed ROI & this will calculate the beta series correlation
% btwn the seed ROI and each voxel in the mask.
% The corr coefficients are Fisher's r-to-z transformed to calculate ttests.
% Will output a ttest nifti file for the stim strings given

% ttest nifti files have 2 volumes:
% 1st vol has the mean r corr coefficient
% & 2nd vol has tstat, computed on Fisher transformed r vals from vol 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% define files etc.

clear all
close all

inDir = ['/home/kelly/ShockAwe/data/results_per_trial'];  % relative to subject directory
outDir = ['/home/kelly/ShockAwe/data/ttests_per_trial'];  % relative to subject directory

maskFile = '/home/kelly/ShockAwe/data/ROIs_tlrc/group_mask.nii.gz';

stims = {'shock','neutral'};

N = 18;

omitProp = zeros(N,length(stims));  % subjects as rows, stim as columns, value denotes proportion of outlier trials excluded

omitOutliers = 1;  % omit outlier beta values? 1 for yes 0 for no
oThresh = 3;

seedStr = 'rVTA_3'; % string to use for out files

% set to 1 to get seed betas using an roi mask, otherwise, set to 0 to get
% seed betas from a text file in the betas directory
useSeedMask = 1; 
seedMaskFile = 'clust_s-n_svc_rVTA.nii.gz'; % give file name if necessary


%% get seed betas

% or to use a seed mask:
if (useSeedMask)
    seedMask = readFileNifti(['/home/kelly/ShockAwe/data/ROIs_tlrc/',seedMaskFile]);
    
else
    
    % to use seed betas from a txt file:
    for c = 1:length(stims)
        seedBs{c}=dlmread(['/home/kelly/ShockAwe/data/betas/',seedStr,'_',stims{c},'_per_trial_betas']);
    end
    
end

%% compute beta series correlations

mask = readFileNifti(maskFile); % get mask file
mask.data = double(mask.data);
idx = find(mask.data);

for c = 1:length(stims)
    
    cd(inDir);
    f = dir(['sa*',stims{c},'*nii*']);
    N = length(f)
    
    for i = 1:N
        
        
        % load betas file
        nii = readFileNifti(f(i).name);
        subData = reshape(nii.data,prod(nii.dim(1:3)),[]); % each volume of voxels is in a row now
         
        
        % get seed betas
        if(useSeedMask)
            
              all_seed_betas = subData(seedMask.data~=0,:); 
           
              w = seedMask.data(seedMask.data~=0); % for weighted mean
              w = repmat(w,1,size(all_seed_betas,2));
          
            if omitOutliers
            
                oidx=find(abs(zscore(all_seed_betas,0,2))>3);
                all_seed_betas(oidx) = nan;
                w(oidx) = nan;
                oCount(i,c) = length(oidx);
            
            end
            
            seedBetas = (nansum(w.*all_seed_betas)./nansum(w))';
                
        else
            
            seedBetas = seedBs{c}(i,~isnan(seedBs{c}(i,:)))';
        
        end
        
        
        % get voxel betas
        voxBetas = subData(mask.data~=0,:)';  % get only voxels within the mask & flip dim
        
        if omitOutliers
            
            gi = logical(abs(zscore(seedBetas))<3); % "good" index - 1 for trials to use and 0 for trials to omit
            oCount2(i,c) = length(find(gi==0));
            sub_rs = corr(seedBetas(gi),voxBetas(gi,:));  % calculate correlation
            
            % recalculate correlation for voxels w/outliers
            [omitTIdx,vIdx]=find(abs(zscore(voxBetas))>3); % get a trial & vox idx for voxel beta outliers
            for j = 1:length(vIdx)
                this_gi = gi; % "good" index for this voxel
                this_gi(omitTIdx(j)) = 0; % omit outlier trial for this voxel
                sub_rs(vIdx(j)) = corr(seedBetas(this_gi),voxBetas(this_gi,vIdx(j)));
            end
            
        else
            
            sub_rs = corr(seedBetas,voxBetas);  % calculate correlation
            
        end % omitOutliers
    
    r{c}(i,:) = sub_rs; % cell array of corr coeffs for each condition
    
end % subjects

z{c} = .5.*log((1+r{c})./(1-r{c}));


%%  ttests
% 
%     meanStat = nanmean(r{c});
%     [~,p,~,stats] = ttest(z{c}); % t-tests on z transformed values
% 
%     meanStatMap = mask.data;
%     meanStatMap(idx) = meanStat;
% 
%     tstatMap = mask.data;
%     tstatMap(idx) = stats.tstat;
% 
%     cd(outDir);
%     outNii = makeGlmNifti(mask,[stims{c},'_',seedStr],meanStatMap,tstatMap);
%     writeFileNifti(outNii);

    
end % stims

cd(outDir);


%% ttests for the contrasts of conditions

% stim 1 - stim 2
meanStat = nanmean(r{1}-r{2});
[~,p,~,stats] = ttest(z{1}-z{2}); % t-tests on z transformed values

meanStatMap = mask.data;
meanStatMap(idx) = meanStat;

tstatMap = mask.data;
tstatMap(idx) = stats.tstat;

outNii = makeGlmNifti(mask,[stims{1},'-',stims{2},'_',seedStr],meanStatMap,tstatMap);
writeFileNifti(outNii);


