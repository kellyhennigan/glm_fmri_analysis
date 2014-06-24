function stats = glm_fmri_fit(Y,X,regIdx)

%%%%%%%%%%%%%%%%%%%%%% INPUTS:

% Y - time series w/n observations in rows (can be 1 or more time series)
% X - design matrix of regressors
% regIdx - (optional) column vector w/length equal to the # of regressors;
% 0 for baseline and >=1 for regressors of interest

% note: # of rows in the time series & model must equal each other

%%%%%%%%%%%%%%%%%%%%%% OUTPUTS:

% a stats struct with the following fields:

% 'B'       - regressor coefficients estimated w/OLS
% 'seB'     - standard error of each beta coefficient
% 'tB'      - this returns the t or F-stat grouped by regIdx as appropriate
% 'pB'      - corresponding p-values for each t-stat
% 'Rsq'     - Rsquared for the full model 
% 'err_ts'   - residual time series


% and if regIdx is provided, also:

% 'Fstat'   - F-stat testing whether the full model is better than baseline
% 'pF'      - p value on F-statistic

% kelly Nov 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define some useful variables 

[n,p] = size(X);  % # of observations and model parameters

% degrees of freedom
df_t = n - 1;     % total degrees of freedom
df_m = p - 1;     % degrees of freedom for model 
df_e = n - p;     % degrees of freedom for error


if size(Y,1)~=n        % if the data and model # of rows aren't equal,
        error('# of rows in data and model must be equal');
end

% if regIdx isnt given, set it to empty
if ~exist('regIdx','var')
    regIdx = [];
end

%% fit the model to the data

fprintf('\nfitting the model...\n');

B = pinv(X'*X)*X'*Y;    % fit model to data with OLS

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

pB = 1 - tcdf(abs(tB), n-2);  % p-value for t-stat

Rsq = 100 * (1 - (SSe ./ SSt) ); % Rsquared; coefficient of multiple determination

fprintf('done.\n');


%% define a stats struct for output with reshaped stats

stats = struct;

stats.B = B;            % beta estimates
stats.seB = seB;        % standard error of beta estimates
stats.tB = tB;          % t-stats for betas
stats.pB = pB;          % p-values for beta t-stats
stats.Rsq = Rsq;        % Rsquared of full model
stats.err_ts = err_ts;  % residuals
stats.df_err = df_e;    % error degrees of freedom


%% if regIdx is given, test the null hypothesis that the full model is 
% no better than the baseline model

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

%  where F is distributed with (? params, df_error) degrees of freedom


if isequal(size(regIdx,2),p)  % if regIdx is given, conduct an F-test
    
    % fit baseline model to data 
    X0 = X(:,regIdx==0);        % baseline model
    B0 = pinv(X0'*X0)*X0'*Y;    
    Yhat0 = X0*B0; 
    SSe0 = sum( (Y - Yhat0).^2 );  % SSe for baseline model
    
    delta_p = count(regIdx~=0); % ? params from baseline > full model 
    Fstat = ((SSe0 - SSe)/delta_p) / MSe; % (? SSe / ? params from reduced > full model) / MSe (full model)
    pF = 1-fcdf(abs(Fstat),delta_p,df_e);  % F is distributed with (? params, df_error) degrees of freedom

    % add some stats to the output structural array 
    stats.Fstat = Fstat;
    stats.pF = pF;

    
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
    
    
end





