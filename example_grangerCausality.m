% Sample neural activity data
% Replace these with your actual data
upper = randn(1, 100);  % Replace with your data
middle = randn(1, 100);  % Replace with your data
deep = randn(1, 100);  % Replace with your data

% Concatenate data for Granger causality analysis
data_for_granger = [upper', middle', deep'];

% Granger causality analysis
max_lag = 5;  % You can adjust this based on your needs
alpha = 0.05; % Significance level

% Function for Granger causality
granger_causality = @(x, y, max_lag, alpha) granger_causality_test(x, y, max_lag, alpha);

% Arrays to store F-statistic values
fStats_upper_middle = zeros(1, 100); % Assuming 100 comparisons
fStats_middle_deep = zeros(1, 100);
fStats_upper_deep = zeros(1, 100);

for i = 1:100  % Perform 100 comparisons, adjust as needed
    % Perform Granger causality analysis between upper and middle
    [pValue, fStat] = granger_causality(data_for_granger(:, 1), data_for_granger(:, 2), max_lag, alpha);
    fStats_upper_middle(i) = fStat;
    
    % Perform Granger causality analysis between middle and deep
    [pValue, fStat] = granger_causality(data_for_granger(:, 2), data_for_granger(:, 3), max_lag, alpha);
    fStats_middle_deep(i) = fStat;
    
    % Perform Granger causality analysis between upper and deep
    [pValue, fStat] = granger_causality(data_for_granger(:, 1), data_for_granger(:, 3), max_lag, alpha);
    fStats_upper_deep(i) = fStat;
end

% Create a bar plot
x = categorical(["Upper vs. Middle", "Middle vs. Deep", "Upper vs. Deep"]);
y = [mean(fStats_upper_middle), mean(fStats_middle_deep), mean(fStats_upper_deep)];  % Use mean for simplicity

bar(x, y);
title('Average Granger Causality F-statistic');
xlabel('Comparison');
ylabel('Average F-statistic');

%%
% Pearson correlation analysis
correlation_upper_middle = corrcoef(upper, middle);
correlation_middle_deep = corrcoef(middle, deep);
correlation_upper_deep = corrcoef(upper, deep);

fprintf('Pearson Correlation:\n');
fprintf('Upper vs. Middle: %.4f\n', correlation_upper_middle(1, 2));
fprintf('Middle vs. Deep: %.4f\n', correlation_middle_deep(1, 2));
fprintf('Upper vs. Deep: %.4f\n', correlation_upper_deep(1, 2));

% Create a bar plot
x = categorical(["Upper vs. Middle", "Middle vs. Deep", "Upper vs. Deep"]);
y = [correlation_upper_middle(1,2), correlation_middle_deep(1,2), correlation_upper_deep(1,2)];  % Use mean for simplicity

bar(x, y);
title('Pearson Correlation');
xlabel('Comparison');
ylabel('Correlation coeffecient');


%% Custom Granger causality function
function [pValue, fStat] = granger_causality_test(x, y, max_lag, alpha)
    % Make sure x & y are the same length
    if (length(x) ~= length(y))
        error('x and y must be the same length');
    end
    
    % Make sure x is a column vector
    [a, b] = size(x);
    if (b > a)
        % x is a row vector -- fix this
        x = x';
    end
    
    % Make sure y is a column vector
    [a, b] = size(y);
    if (b > a)
        % y is a row vector -- fix this
        y = y';
    end
    
    % Make sure max_lag is >= 1
    if max_lag < 1
        error('max_lag must be greater than or equal to one');
    end
    
    % First find the proper model specification using the Bayesian Information
    % Criterion for the number of lags of x
    T = length(x);
    BIC = zeros(max_lag, 1);
    
    % Specify a matrix for the restricted RSS
    RSS_R = zeros(max_lag, 1);
    i = 1;
    while i <= max_lag
        ystar = x(i + 1:T, :);
        xstar = [ones(T - i, 1) zeros(T - i, i)];
        
        % Populate the xstar matrix with the corresponding vectors of lags
        j = 1;
        while j <= i
            xstar(:, j + 1) = x(i + 1 - j:T - j);
            j = j + 1;
        end
        
        % Apply the regress function. b = betahat, bint corresponds to the 95%
        % confidence intervals for the regression coefficients and r = residuals
        [b, ~, r] = regress(ystar, xstar);
        
        % Find the Bayesian information criterion
        BIC(i, :) = T * log(r' * r / T) + (i + 1) * log(T);
        
        % Put the restricted residual sum of squares in the RSS_R vector
        RSS_R(i, :) = r' * r;
        
        i = i + 1;
    end
    
    [~, x_lag] = min(BIC);
    
    % First find the proper model specification using the Bayesian Information
    % Criterion for the number of lags of y
    BIC = zeros(max_lag, 1);
    
    % Specify a matrix for the unrestricted RSS
    RSS_U = zeros(max_lag, 1);
    i = 1;
    while i <= max_lag
        ystar = x(i + x_lag + 1:T, :);
        xstar = [ones(T - (i + x_lag), 1) zeros(T - (i + x_lag), x_lag + i)];
        
        % Populate the xstar matrix with the corresponding vectors of lags of x
        j = 1;
        while j <= x_lag
            xstar(:, j + 1) = x(i + x_lag + 1 - j:T - j, :);
            j = j + 1;
        end
        
        % Populate the xstar matrix with the corresponding vectors of lags of y
        j = 1;
        while j <= i
            xstar(:, x_lag + j + 1) = y(i + x_lag + 1 - j:T - j, :);
            j = j + 1;
        end
        
        % Apply the regress function. b = betahat, bint corresponds to the 95%
        % confidence intervals for the regression coefficients and r = residuals
        [b, ~, r] = regress(ystar, xstar);
        
        % Find the Bayesian information criterion
        BIC(i, :) = T * log(r' * r / T) + (i + 1) * log(T);
        
        RSS_U(i, :) = r' * r;
        
        i = i + 1;
    end
    
    [~, y_lag] = min(BIC);
    
    % The numerator of the F-statistic
    F_num = ((RSS_R(x_lag, :) - RSS_U(y_lag, :)) / y_lag);
    
    % The denominator of the F-statistic
    F_den = RSS_U(y_lag, :) / (T - (x_lag + y_lag + 1));
    
    % The F-Statistic
    F = F_num / F_den;
    
    % Critical value from F-distribution
    c_v = finv(1 - alpha, y_lag, (T - (x_lag + y_lag + 1)));
    
    % Compute p-value
    pValue = 1 - fcdf(F, y_lag, (T - (x_lag + y_lag + 1)));
    
    % Output F-statistic and p-value
    fStat = F
end