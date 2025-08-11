function W_score=parcelwise_Wscore(HC_data,Case_data,Covar_hc,Covar_case)

% check the dimensions
[~, n_Parcel] = size(HC_data);
[n_Cases, ~] = size(Case_data);

%% W‐score maps
% Create an  array to save the w-score result  
W_score = zeros(n_Cases, n_Parcel); 

% load cortical thickness data of HCs
for p=1:n_Parcel
    yHC = HC_data(:, p);
    %Build a model based on data from HCs
    lm = fitlm(Covar_hc, yHC, 'linear');
    %Perform predictions for the HCs
    yHC_pred = predict(lm, Covar_hc);
    %Calculate residuals for the HC group
    residuals_hc = yHC - yHC_pred;
    %Calculate the standard deviation of HC residuals
    sd_hc = std(residuals_hc, 0, 1);
    if sd_hc < 1e-12
        sd_hc = 1e-12;  
    end
    %Predict CT for the cases using the model trained on HC
    yCases_pred = predict(lm, Covar_case);%n_cases x 1
    %True values for the cases 
    yCases = Case_data(:, p);
    %Calculate w-score map, W‐score = (actual – expected)/SD
    W_score(:, p) = (yCases - yCases_pred) ./ sd_hc;
end