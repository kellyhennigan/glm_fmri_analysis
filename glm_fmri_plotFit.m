function fig = glm_fmri_plotFit(ts,X,regIdx)

% this function plots time series data and a model fit to the data. 

%%%%% INPUTS

% ts      - raw fmri time series (from a voxel or an roi)
% X       - glm model; # of rows should match ts; 1 column for each regressor
% regIndx (optional) - index specifying regressors of interest (>0) and nuisance
%                      regressors (0)

%%%%% OUTPUTS 

% a figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cb=cbrewer('qual', 'Set1', 8);
% colormap(cb)

% set regIdx to be a vector of zeros if not defined
if notDefined('regIdx')
    regIdx=zeros(1,size(X,2));
end


colors = solarizedColors(size(X,2));

B = glm_fmri_fit(ts,X,regIdx,'B');

Yhat = X*B;

% B = pinv(X'*X)*X'*ts;

TRs = 1:length(ts);


% j=sum(X(:,regIdx==1),2).*mean(B(regIdx==1)); % juice regressor * beta
% n=sum(X(:,regIdx==2),2).*mean(B(regIdx==2)); % neutral ""
% s=sum(X(:,regIdx==3),2).*mean(B(regIdx==3)); % shock ""


fig = figure
set(gcf,'Position',[1 300 1360 384])
plot(TRs,ts','k')
hold on
plot(TRs,Yhat','b') % total model fit
plot(TRs,j+mean(ts),'color',colors(1,:)) % juice
plot(TRs,n+mean(ts),'color',colors(2,:)) % neutral
plot(TRs,s+mean(ts),'color',colors(3,:)) % shock
legend('time series','model fit','juice','neutral','shock')
% subplot(2,1,2)
% set(gcf,'DefaultAxesColorOrder',colors)
% h=plot(TRs,X(:,regIdx~=0)')
% legend('juice','neutral','shock')






