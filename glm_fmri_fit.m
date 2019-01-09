function stats = glm_fmri_fit(Y,X,regIdx,statStr)

%%%%%%%%%%%%%%%%%%%%%% INPUTS:

% Y - time series w/n observations in rows (can be 1 or more time series)

% X - design matrix of regressors. ***AT LEAST 1 COLUMN MUST CONTAIN ALL
% CONSTANTS TO MODEL MEAN ACTIVITY***


% regIdx - (optional) column vector w/length equal to the # of regressors;
%          0 for baseline and >=1 for regressors of interest

% statStr - string specifying a stat to return (e.g., 'B', or 'Rsq'). If
%           not given, then a structural array containing all stats will be
%           returned


% note: # of rows in the time series & model must equal each other

%%%%%%%%%%%%%%%%%%%%%% OUTPUTS:

% either a specific stat as specified with input 'statStr', or a structural
% array that contains all computed stats as fields:

% 'B'       - regressor coefficients estimated w/OLS
% 'seB'     - standard error of each beta coefficient
% 'tB'      - this returns the t or F-stat grouped by regIdx as appropriate
% 'pB'      - corresponding p-values for each t-stat
% 'Rsq'     - Rsquared for the full model (ordinary R2)
% 'err_ts'   - residual time series

% and if regIdx is provided, also:

% 'Fstat'   - F-stat testing whether the full model is better than baseline
% 'pF'      - p value on F-statistic

%%%%%%%%%%%%%%%%%%%%%% EXAMPLES:

% stats = glm_fmri_fit(Y,X);
% stats = glm_fmri_fit(Y,X,regIdx);
% pF = glm_fmri_fit(Y,X,regIdx,'pF');
% tB = glm_fmri_fit(Y,X,[],'tB');

% kelly Nov 2013

% see: http://www.fil.ion.ucl.ac.uk/spm/doc/books/hbf2/pdfs/Ch8.pdf
%   for more info on fmri model fitting

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some useful variables

[n,p] = size(X);  % N observations and p model parameters

if size(Y,1)~=n        % if the data and model # of rows aren't equal,
    error('# of rows in data and model must be equal');
end

% if regIdx isnt given, set it to empty
if ~exist('regIdx','var')
    regIdx = [];
end

if ~exist('statStr','var')
    statStr = '';
end


%% fit the model to the data

% degrees of freedom
df_t = n - 1;     % total degrees of freedom
df_m = p - 1;     % degrees of freedom for model
df_e = n - p;     % degrees of freedom for error


fprintf('\nfitting the model...\n');

B = pinv(X'*X)*X'*Y;    % pinv(x) gives the (Moore-Penrose) pseudo inverse of x

Yhat = X*B;             % model's prediction for Y

err_ts = Y-Yhat;        % error time series

SSt = sum( (Y - repmat(mean(Y),n,1)).^2 );   % sum of squares (total)

SSm = sum( (Yhat - repmat(mean(Y),n,1)).^2 );      % sum of squares (model)

SSe = sum( (Y - Yhat).^2 );   % sum of squares (error)

MSt = SSt./df_t;              % mean sq error aka sample variance

MSm = SSm./df_m;              % mean sq error aka explained variance

MSe = SSe./df_e;              % mean sq error aka unexplained variance

seB = sqrt(diag(pinv(X'*X))*MSe);   % standard error of beta coefficients

tB = B ./ seB;     % distributed as t w/ n-2 df

pB = (1 - tcdf(abs(tB), df_e))./2;  % p-value for t-stat

Rsq = 1 - (SSe ./ SSt) ; % Rsquared; coefficient of multiple determination

fprintf('done.\n');


%% if regIdx is given, test the null hypothesis that the full model
% is no better than the baseline model

% full model = baseline + regressors of interest
% baseline model         = X(regIdx==0)
% regressors of interest = X(regIdx~=0)

%          H0:   B1 = B2 = ... = Bj = 0; (for regressors of interest)
%          H1:   Bj ~= 0 for at least 1 regressor of interest

% Appropriate stat for this test is the F-statistic.

% To test the full model against a constant model (e.g, the situation where
% the baseline model is just one constant regressor):

%              MSm         explained variance
%      F  =   -----   =   --------------------
%              MSe        unexplained variance

%  where F is distributed with (df_model, df_error) degrees of freedom



% To test a full model that has p total regressors against a reduced model
% of only baseline regressors:

%             ? SSe / ? params from reduced > full model
%      F  =   -----------------------------------------
%                        MSe (full model)

%  where F is distributed with (p -1 , n - p) degrees of freedom

% or equivalently:

%                           Rsq * (n - p)
%   F(n-1,n-p)  =   -----------------------------------------
%                        (1 - Rsq) * (p - 1)

if isequal(size(regIdx,2),p)  % if regIdx is given, conduct an F-test
    
    % fit baseline model to data
    X0 = X(:,regIdx==0);        % baseline model
    B0 = pinv(X0'*X0)*X0'*Y;
    Yhat0 = X0*B0;
    SSe0 = sum( (Y - Yhat0).^2 );  % SSe for baseline model
    
    delta_p = sum(regIdx~=0); % # of params in full - baseline model
    Fstat = ((SSe0 - SSe)./repmat(delta_p,1,size(SSe,2))) ./ MSe; % (? SSe / ? params from reduced > full model) / MSe (full model)
    pF = 1-fcdf(abs(Fstat),delta_p,df_e);  % F is distributed with (? params, df_error) degrees of freedom
    
    %
    %     df_top = length(find(regIdx~=0)); % df_top is for numerator; df is for denominator
    %
    %     Yhat0 = X(:,regIdx==0)*B(regIdx==0,:); % baseline model's prediction of data
    %
    %     SSe0 = sum( (Y - Yhat0).^2 );        % sum of squares error for baseline model
    %
    %     Fstat = (SSe0 - SSe)./df_top ./ MSe; % F statistic
    %
    %     pF = 1-fcdf(abs(Fstat),df_top,df);
    
    
else
    
    % set F-stat and corresponding p value to nan if regIdx isn't given
    Fstat = nan;
    pF = nan;
    
end


%% which stats to return as output?

switch statStr
    
    case 'B'
        stats = B;
        
    case 'seB'
        stats = seB;
        
    case 'tB'
        stats = tB;
        
    case 'pB'
        stats = pB;
        
    case 'Rsq'
        stats = Rsq;
        
    case 'err_ts'
        stats = err_ts;
        
    case 'df_err'
        stats = df_e;
        
    case 'Fstat'
        stats = Fstat;
        
    case 'pF'
        stats = pF;
        
    otherwise
        
        stats = struct;
        
        stats.B = B;            % beta estimates
        stats.seB = seB;        % standard error of beta estimates
        stats.tB = tB;          % t-stats for betas
        stats.pB = pB;          % p-values for beta t-stats
        stats.Rsq = Rsq;        % Rsquared of full model
        stats.err_ts = err_ts;  % residuals
        stats.df_err = df_e;    % error degrees of freedom
        
        stats.Fstat = Fstat;   % F stat testing full vs baseline model (as defined by regIdx)
        stats.pF = pF;         % p value for F stat
        
end % switch statStr

end % function





