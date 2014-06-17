function regTS = glm_fmri_createRegTS(onsets,TR,nt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this function creates regressor time series for fmri analysis.


% INPUTS:
%   onsets - M x N matrix of onset times (in seconds) of events of
%           interest; time should be in seconds and relative to the scan
%           start time. Onset times within each column will be modeled
%           as a single regressor time series. Negative #s, zero, and nan
%           values are ignored. 

%   TR - scan repetition time (in seconds)
%   nt - # of volumes acquired in the scan run being modeled.
           
% OUTPUTS:
%     regTS - nt x N regressor times series where T is the number of modeled
%           time points and N is the number of regressors.


% NOTES: 


% kjh, 20-Nov-2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% (upsampled) sampling rate for convolution purposes, in seconds. Probably shouldn't be changed
dt = .1;  
hrf = spm_hrf(dt); % hemodynamic response function w/spm's default parameters
hrf = hrf./max(hrf); % scale it so that the max value is 1


% change the unit of time to .1 sec 
onsets = round(onsets./dt);


% define regressor time series as vectors of zeros w/ ones at event onsets
ts = zeros(nt*TR/dt,size(onsets,2));

for ic = 1:size(onsets,2)  % for each column of onset times
    ts(onsets(:,ic),ic) = ones; % set onsets in time series to ones
    reg_ts(:,ic) = conv(ts(:,ic),hrf); % convolve regressor time series with the hrf
end

% downsample the time series to match the fmri data sampling rate (TR)
regTS = reg_ts(1:TR/dt:nt(1)*TR/dt,:);


end
