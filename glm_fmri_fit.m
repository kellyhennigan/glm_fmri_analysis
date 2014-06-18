function stats = glm_fmri_fit(Y,X,regIdx)

%%%%%%%%%%%%%%%%%%%%%% INPUTS:

% Y - n time series of m time points 
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

[nt,np] = size(X);      % # of time points and model parameters
df = nt - np;       % degrees of freedom

if size(Y,1)~=nt        % if the data and model # of rows aren't equal,
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

SSt = sum( (Y - repmat(mean(Y),nt,1)).^2 );   % sum of squares (total)

SSe = sum( (Y - Yhat).^2 );      % sum of squares (error)

MSe = SSe./df;                   % mean sq error aka error variance

seB = sqrt(diag(pinv(X'*X))*MSe);   % standard error of beta coefficients

tB = B ./ seB;     % distributed as t w/ nt-2 df

pB = 1 - tcdf(abs(tB), nt-2);  % p-value for t-stat

Rsq = 100 * (1 - (SSe ./ SSt) ); % Rsquared; coefficient of multiple determination

fprintf('done.\n');


%% define a stats struct for output with reshaped stats

stats = struct;

stats.B = B; 
stats.seB = seB;
stats.tB = tB;
stats.pB = pB;
stats.Rsq = Rsq;
stats.err_ts = err_ts;
stats.df = df;

%% if regIdx is given, test the null hypothesis that the full model is 
% no better than the baseline model

if isequal(size(regIdx),[1,size(X,2)])
    
    df_top = length(find(regIdx~=0)); % df_top is for numerator; df is for denominator
    
    Yhat0 = X(:,regIdx==0)*B(regIdx==0,:); % baseline model's prediction of data
    
    SSe0 = sum( (Y - Yhat0).^2 );        % sum of squares error for baseline model
    
    Fstat = (SSe0 - SSe)./df_top ./ MSe; % F statistic
    
    pF = 1-fcdf(abs(Fstat),df_top,df);
    
    % add some stats to the output structural array 
    stats.Fstat = Fstat;
    stats.pF = pF;
    
end





