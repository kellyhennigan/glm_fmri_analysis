function plotDesignMat(X,regLabels,regIdx)
% --------------------------------
% function to plot a design matrix made using glm_fmri_mat function. Takes
% in a design matrix and separately plots the model of interest and
% baseline model, as specified from regIdx.

% INPUT:
% X - glm design matrix with a column for each regressor
% regLabels - cell array of string labels corresponding to each reg column
% regIdx - vector w/ 1,2, etc. for each stim specified by stims and 0 for
%       all other (baseline) regs
%
%
% OUTPUT:
%   figures plotting the model of interest and baseline model
%

% author: Kelly, kelhennigan@gmail.com, 26-Nov-2014

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xoi = X(:,regIdx~=0);    % model of interest

figure
imagesc(Xoi, [0 1])
title('GLM Design Matrix - regressors of interest')
set(gca,'XTick',1:size(Xoi,2))
set(gca,'XTickLabel',regLabels(regIdx~=0))
rotateXLabels(gca,45)
colormap(gray)


Xbase = X(:,regIdx==0);  % baseline model

figure
imagesc(Xbase, [0 1])
title('GLM Design Matrix - baseline model')
set(gca,'XTick',1:size(Xbase,2))
set(gca,'XTickLabel',regLabels(regIdx==0))
rotateXLabels(gca,45)
colormap(gray)


