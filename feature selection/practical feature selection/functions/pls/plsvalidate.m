function [Yp, RMSE] = plsvalidate(model, varargin)

%DESCRIPTION:
%This function validates PLS models previously trained using 'plstrain.m'. 
%If the number of LVs were optimized for that model (and not manually set), 
%this function re-optimizes the LVs for each validation and thus performs 
%a double cross-validation. The following names are given to the validation 
%layers:
%- Calibration: Inner most set used to train the models (for each number of LVs).
%- Validation: Middle set used to test the models for each number of LVs.
%- Prediction: Outer most set used to test the model with the number of LVs
%  that gave the highest performance on the validation set.
%The function can also used to perform a permutation test by permuting the 
%samples in each iteration. For more information on double (cross)-validation, 
%see Smit et al, 2007, 'Assessing the statistical validity of proteomics based 
%biomarkers'.
%
%INPUT:
%- model = The PLS model to validate, obtained using 'plstrain.m'.
%
%- The following name-value pairs can be used:
%  - 'scheme': Validation scheme to use for the validation:
%    - 'random' = Validation using a single random subset.
%    - 'kfold' = Random cross-validation using k folds.
%    - 'vblindsi' = Venetian blinds cross-validation on the sampling index
%      (recommended for time-series data).
%    - 'vblindsy'= Venetian blinds cross-validation on the first column of Y 
%      (default).
%    - 'leaveout' Leave-one-out cross-validation.
%    - 'cblock' = Contiguous block cross-validation.
%  - 'n': Validation repeats to use (for non-deterministic schemes, 
%     default = 1).
%  - 'k': For cross-validation, this is the number of partitions the data is
%    divided into (default = 5). For single test set validation, a fraction of 
%    1/k is taken out of the data for validation (so for k=5, 20% of the data 
%    is used for testing.
%  - 'alpha': The confidence level that is used to calculate the two-sided
%    confidence interval of the prediction in case the validation is repeated. 
%    An alpha of 0.05 corresponds to a 95%-confidence interval. The output
%    is given as 'Ypci', which respecticely gives the mean, standard
%    deviation and lower and upper confidence limits.
%  - 'perm': Whether to perform a permutation test.
%  - 'layer': Force the function to perform a 'single' or 'double'
%    validation, regardless of how the actual model is optimized.
%  - 'optscheme': Validation scheme used for the optimization of the model (in 
%    case of double validation).
%  - 'optn': Validation repeats used for the optimization of the model (in 
%    case of double validation).
%  - 'optk': Validation folds used for the optimization of the model (in 
%    case of double validation).
%  - 'optlv': Maximum number of LVs to consider for optimization of the
%    model (in case of double validation).
%
%OUTPUT:
%- Yp: The validated Y-value, corresponding to the outermost validation
%  layer (Yval in case of single validation, Ypred in case of double
%  validation.
%- RMSE: The root mean squared error between prediction and reference for
%  the outermost validation layer. For discrimination and classification
%  models, this is calculated before transforming the continuous PLS output
%  to the dummy matrix.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 2021
%
%SYNTAX:
%[Yp, RMSE] = plsvalidate(model, ...)
%[validation] = plsvalidate(model, ...)
%[validation] = plsvalidate(model, 'scheme', 'kfold', 'n', 10, 'k', 5, ...)

%Parse input:
scheme = [];
n = [];
k = [];
alpha = [];
perm = false;
layer = [];
optscheme = [];
optn = [];
optk = [];
optlv = [];
for i=1:2:length(varargin)
    eval([varargin{i} ' = varargin{i+1};']);
end

%Correct input:
if k<1
    k = ceil(1/k);
end
if optk<1
    optk = ceil(1/optk);
end

%Find settings for inner model optimization (if applicable):
if isempty(layer)
    if isfield(model, 'optimization')
        layer = 'double';
    else
        layer = 'single';
    end
else
    layer = 'single';
end

%Find (and copy) settings in case of single cross-validation
if strcmp(layer, 'single')
    
    %If no scheme is given, use the default:
    if isempty(scheme)
        scheme = 'vblinds';
    end
    if isempty(k)
        k = 5;
    end
    if isempty(n)
        n = 1;
    end
end
    
%find (and copy) settings in case of double cross-validation
if strcmp(layer, 'double')
    
    %If no outer scheme is given, copy the one used for model optimization
    if isempty(scheme)
        scheme = model.optimization.optscheme;
    end
    if isempty(k)
        k = model.optimization.optk;
    end
    if isempty(n)
        n = model.optimization.optn;
    end
    
    %If no inner scheme is give, copy the outer scheme
    if isempty(optscheme)
        optscheme = scheme;
    end
   if isempty(optk)
        optk = k;
    end
    if isempty(optn)
        optn = n;
    end
    if isempty(optlv)
        if isfield(model, 'optimization')
        	optlv = model.optimization.LVmax;
        else
        	optlv = model.LV;
        end
    end
end

%Determing validation set(s):
if strcmp(scheme, 'vblindsi') || strcmp(scheme, 'vblinds')
    n = 1;
    test = false(n, k, size(model.Y, 1));
    for i=1:size(test, 2)
        test(1, i, i:size(test, 2):size(test, 3)) = true;
    end
elseif strcmp(scheme, 'vblindsy')
    n = 1;
    test = false(n, k, size(model.Y, 1));
    [~, s] = sort(model.Y(:, 1));
    for i=1:size(test, 2)
        test(1, i, s(i:size(test, 2):size(test, 3))) = true;
    end
elseif strcmp(scheme, 'kfold')
    test = false(n, k, size(model.Y, 1));
    for in = 1:size(test, 1)
        s = randperm(size(test, 3));
        for i=1:size(test, 2)
            test(in, i, s(i:size(test, 2):size(test, 3))) = true;
        end
    end
elseif strcmp(scheme, 'random')
    test = false(n, 1, size(model.Y, 1));
    for in = 1:size(test, 1)
        s = randperm(size(test, 3));
        test(in, 1, s(1:floor(size(test, 3) * (1/k)))) = true;
    end
elseif strcmp(scheme, 'leaveout')
    n = 1;
    test = false(n, k, size(model.Y, 1));
    k = size(model.Y, 1);
    test = false(1, size(model.Y, 1), size(model.Y, 1));
    test(1, find(eye(size(model.Y, 1)))) = true;
elseif strcmp(scheme, 'cblock')
    n = 1;
    test = false(n, k, size(model.Y, 1));
	for i=1:size(test, 2)
        test(n, i, floor((i-1) * (size(test, 3) / size(test, 2)) + 1):floor((i) * (size(test, 3) / size(test, 2)))) = true;
    end
else
    error('Invalid validation scheme. Please check the help-information and your input.');
end

%Initialize validation results:
Ypcal = NaN([size(test, 1), size(test, 2), size(model.Y, 1), size(model.Y, 2)]);
Ypval = NaN([size(test, 1), size(test, 2), size(model.Y, 1), size(model.Y, 2)]);
Yppred = NaN([size(test, 1), size(test, 2), size(model.Y, 1), size(model.Y, 2)]);

%For each validation repeat:
for in = 1:size(test, 1)
    
    %Make a temporary copy for the model:
    modeltemp = model;
    
    %Permute the data, if asked to do so:
    if perm==true
        modeltemp.X = modeltemp.X(randperm(size(modeltemp.X, 1)), :);
    end
    
    %Model for each fold:
    for ik = 1:size(test, 2)
        
        %Split data into training and test set:
        Xtrain = modeltemp.X(~test(in, ik, :), :);
        Ytrain = modeltemp.Y(~test(in, ik, :), :);
        Xtest = modeltemp.X(test(in, ik, :), :);
        Ytest = modeltemp.Y(test(in, ik, :), :);
        
        %Perform single validation:
        if strcmp(layer, 'single')
            temp = plstrain(Xtrain, Ytrain, modeltemp.Xpp, modeltemp.Ypp, modeltemp.LV);
            [Ypcal(in, ik, ~test(in, ik, :), :), ~] = plspredict(temp, Xtrain);
            [Ypval(in, ik, test(in, ik, :), :), ~] = plspredict(temp, Xtest, [], 'regression');
        
        %Perform double validation:
        elseif strcmp(layer, 'double')
            temp = plstrain(Xtrain, Ytrain, modeltemp.Xpp, modeltemp.Ypp, -abs(optlv), 'optscheme', optscheme, 'optn', optn, 'optk', optk);
            Ypcal(in, ik, ~test(in, ik, :), :) = temp.optimization.Ypcal;
            Ypval(in, ik, ~test(in, ik, :), :) = temp.optimization.Ypval;
            [Yppred(in, ik, test(in, ik, :), :), ~] = plspredict(temp, Xtest, [], 'regression');
        end
    end
end

%Average results over validation folds:
if size(Ypcal, 1)==1
    Ypcal = reshape(squeeze(nanmean(Ypcal, 2)), size(model.Y));
    A = NaN([1, size(Ypcal)]);
    A(1, :, :) = Ypcal;
    Ypcal = A;
else
    Ypcal = reshape(squeeze(nanmean(Ypcal, 2)), [n size(model.Y)]);
end
if size(Ypval, 1)==1
    Ypval = reshape(squeeze(nanmean(Ypval, 2)), size(model.Y));
    A = NaN([1, size(Ypval)]);
    A(1, :, :) = Ypval;
    Ypval = A;
else
    Ypval = reshape(squeeze(nanmean(Ypval, 2)), [n size(model.Y)]);
end
if size(Yppred, 1)==1
    Yppred = reshape(squeeze(nanmean(Yppred, 2)), size(model.Y));
    A = NaN([1, size(Yppred)]);
    A(1, :, :) = Yppred;
    Yppred = A;
else
    Yppred = reshape(squeeze(nanmean(Yppred, 2)), [n size(model.Y)]);
end

%Calculate confidence intervals of results:
if ~isempty(alpha) && size(test, 1)>1

    %If alpha is higher than 1, the user probably inputed a CI-percentage that 
    %needs to be converted to an alpha value:
    if alpha>1
        alpha = 1-alpha./100;
    end
    
    %Look up Z-value for the alpha-value:
    ztable = [0,Inf;0.01,2.576;0.02,2.326;0.03,2.170;0.04,2.054;0.05,1.960;0.06,1.881;0.07,1.812;0.08,1.751;0.09,1.695;0.1,1.645;0.11,1.598;0.12,1.555;0.13,1.514;0.14,1.476;0.15,1.440;0.16,1.405;0.17,1.372;0.18,1.341;0.19,1.311;0.2,1.282;0.21,1.254;0.22,1.227;0.23,1.2;0.24,1.175;0.25,1.150;0.26,1.126;0.27,1.103;0.28,1.080;0.29,1.058;0.3,1.036;0.31,1.015;0.32,0.9940;0.33,0.9740;0.34,0.9540;0.35,0.9350;0.36,0.9150;0.37,0.8960;0.38,0.8780;0.39,0.86;0.4,0.8420;0.41,0.8240;0.42,0.8060;0.43,0.7890;0.44,0.7720;0.45,0.7550;0.46,0.7390;0.47,0.7220;0.48,0.7060;0.49,0.69;0.5,0.6740;0.51,0.6590;0.52,0.6430;0.53,0.6280;0.54,0.6130;0.55,0.5980;0.56,0.5830;0.57,0.5680;0.58,0.5530;0.59,0.5390;0.6,0.5240;0.61,0.51;0.62,0.4960;0.63,0.4820;0.64,0.4680;0.65,0.4540;0.66,0.44;0.67,0.4260;0.68,0.4120;0.69,0.3990;0.7,0.3850;0.71,0.3720;0.72,0.3580;0.73,0.3450;0.74,0.3320;0.75,0.3190;0.76,0.3050;0.77,0.2920;0.78,0.2790;0.79,0.2660;0.8,0.2530;0.81,0.24;0.82,0.2280;0.83,0.2150;0.84,0.2020;0.85,0.1890;0.86,0.1760;0.87,0.1640;0.88,0.1510;0.89,0.1380;0.9,0.1260;0.91,0.1130;0.92,0.1;0.93,0.088;0.94,0.075;0.95,0.063;0.96,0.05;0.97,0.038;0.98,0.025;0.99,0.013;1,0];
    alpha = round(alpha, 2);
    z = ztable(find(ztable(:, 1)==alpha, 1), 2);
       
    %Calculate the mean, standard deviation and confidence interval for Ypcal:
    D = size(Ypcal);
    Ypcal_ci = zeros([4, D(2:end)]);
    Ypcal_ci(1, :) = reshape(squeeze(mean(Ypcal, 1)), 1, prod(D(2:end)));
    Ypcal_ci(2, :) = reshape(squeeze(std(Ypcal, 0, 1)), 1, prod(D(2:end)));
    Ypcal_ci(3, :) = reshape(squeeze(Ypcal_ci(1, :)) - (z.*(squeeze(Ypcal_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    Ypcal_ci(4, :) = reshape(squeeze(Ypcal_ci(1, :)) + (z.*(squeeze(Ypcal_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    
    %Calculate the mean, standard deviation and confidence interval for Ypval:
    D = size(Ypval);
    Ypval_ci = zeros([4, D(2:end)]);
    Ypval_ci(1, :) = reshape(squeeze(mean(Ypval, 1)), 1, prod(D(2:end)));
    Ypval_ci(2, :) = reshape(squeeze(std(Ypval, 0, 1)), 1, prod(D(2:end)));
    Ypval_ci(3, :) = reshape(squeeze(Ypval_ci(1, :)) - (z.*(squeeze(Ypval_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    Ypval_ci(4, :) = reshape(squeeze(Ypval_ci(1, :)) + (z.*(squeeze(Ypval_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    
    %Calculate the mean, standard deviation and confidence interval for Yppred:
    D = size(Yppred);
    Yppred_ci = zeros([4, D(2:end)]);
    Yppred_ci(1, :) = reshape(squeeze(mean(Yppred, 1)), 1, prod(D(2:end)));
    Yppred_ci(2, :) = reshape(squeeze(std(Yppred, 0, 1)), 1, prod(D(2:end)));
    Yppred_ci(3, :) = reshape(squeeze(Yppred_ci(1, :)) - (z.*(squeeze(Yppred_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    Yppred_ci(4, :) = reshape(squeeze(Yppred_ci(1, :)) + (z.*(squeeze(Yppred_ci(2, :))./(n^0.5))), prod(D(2:end)), 1);
    
    %Calculate the mean, standard deviation and confidence interval for
    %RMSE (outer most validation):
    RMSEs = zeros(D(1), 1);
    for i=1:D(1)
        if strcmp(layer, 'double')
            RMSEs(i) = nanmean((squeeze(Yppred(i, :))' - model.Y(:)).^2).^0.5;
        else
            RMSEs(i) = nanmean((squeeze(Ypval(i, :))' - model.Y(:)).^2).^0.5;
        end
    end
    RMSE_ci = [mean(RMSEs); std(RMSEs); mean(RMSEs) - z*std(RMSEs)/D(1).^0.5; mean(RMSEs) + z*std(RMSEs)/D(1).^0.5];
    
    %Calculate the mean, standard deviation and confidence interval for the R (in case of regression, and for outer most validation):
    if strcmp(model.type, 'regression')
        Rs = zeros(D(1), 1);
        for i=1:D(1)
            if strcmp(layer, 'double')
                Rs(i) = min(min(corrcoef(Yppred(i, :), model.Y(:))));
            else
                Rs(i) = min(min(corrcoef(Ypval(i, :), model.Y(:))));
            end
        end
        R_ci = [mean(Rs); std(Rs); mean(Rs) - z*std(Rs)/D(1).^0.5; mean(Rs) + z*std(Rs)/D(1).^0.5];
    end
end

%Average results over validation repeats:
Ypcal = reshape(squeeze(nanmean(Ypcal, 1)), size(model.Y));
Ypval = reshape(squeeze(nanmean(Ypval, 1)), size(model.Y));
Yppred = reshape(squeeze(nanmean(Yppred, 1)), size(model.Y));

%Select correct output for calculating (outermost) validated performance:
if strcmp(layer, 'double')
    Yp = Yppred;
else
    Yp = Ypval;
end

%Calculate validated RMSE (for all model types):
RMSE = sqrt(mean((model.Y(:)-Yp(:)).^2));

%Calculate R for regression models:
if strcmp(model.type, 'regression')
    R = min(min(corrcoef(model.Y(:), Yp(:))));
end

%Transform Yp to dummy matrix in case of discrimination or classification:
if strcmp(model.type, 'discrimination')
    Ypcal = double(repmat(max(Ypcal, [], 2), 1, size(Ypcal, 2)) == Ypcal);
    Ypval = double(repmat(max(Ypval, [], 2), 1, size(Ypval, 2)) == Ypval);
    Yppred = double(repmat(max(Yppred, [], 2), 1, size(Yppred, 2)) == Yppred);
    Yp = double(repmat(max(Yp, [], 2), 1, size(Yp, 2)) == Yp);
elseif strcmp(model.type, 'classification')
    Ypcal = double(Ypcal>=0.5);
    Ypval = double(Ypval>=0.5);
    YpYppredcal = double(Ypcal>=0.5);
    Yp = double(Yp>=0.5);
end

%Calculate CM for discrimination and classification models (true classes in rows and predicted classes in column), and PA and MCC if possible:
if strcmp(model.type, 'discrimination') || strcmp(model.type, 'classification')
    CM = zeros(size(model.Y, 2));
    for i=1:size(CM, 1)
        CM(i, :) = sum(repmat(model.Y(:, i), 1, size(model.Y, 2)) .* Yp);
    end
    if (size(model.Y, 2)==1) && (~any(sum(model.Y, 2)~=1))
        TP = CM(1, 1);
        TN = CM(2, 2);
        FP = CM(1, 2);
        FN = CM(2, 1);
        PA = (TP+TN) / (TP+TN+FP+FN);
        if ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) == 0
            MCC = ((TP * TN) - (FP * FN)) / 1;
        else
            MCC = ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));
        end
    end
end

%Prepare validation-structure, if requested to do so:
if nargout==1
	validation = [];
    if strcmp(layer, 'double')
        validation.Yppred = Yppred;
    end
    validation.Ypval = Ypval;
    validation.Ypcal = Ypcal;
    validation.Yref = model.Y;
    validation.X = model.X;
    validation.RMSE = RMSE;
    if exist('R')==1
        validation.R = R;
    end
    if exist('CM')==1
        validation.CM = CM;
    end
    if exist('PA')==1
        validation.PA = PA;
    end
    if exist('MCC')==1
        validation.MCC = MCC;
    end
    validation.layer = layer;
    validation.scheme = scheme;
    validation.n = n;
    validation.k = k;
    validation.permutation = perm;
    if ~isempty(alpha)
        validation.confidence.alpha = alpha;
        if strcmp(layer, 'double')
            validation.confidence.Yppred_ci = Yppred_ci;
        end
        validation.confidence.Ypval_ci = Ypval_ci;
        validation.confidence.Ypcal_ci = Ypcal_ci;
        validation.confidence.RMSE_ci=RMSE_ci;
        if exist('R_ci')==1
            validation.confidence.R_ci = R_ci;
        end
    end
    if strcmp(layer, 'double')
        validation.optimization.valscheme = optscheme;
        validation.optimization.valn = optn;
        validation.optimization.valk = optk;
        validation.optimization.LVmax = optlv;
    end
    validation.info = {
        'Yppred', 'Predicted reference data of prediction (outer/double validation)';
        'Ypval', 'Predicted reference data of validation (inner/single validation)';
        'Ypcal', 'Predicted reference data of calibration';
        'Yref', 'Reference dependent data';
        'X', 'Inputted independent data';
        'RMSE', 'The root mean squared error between prediction and reference for the outermost validation layer. For discrimination and classification models, this is calculated before transforming the continuous PLS output to the dummy matrix';
        'R','Correlation between prediction and reference Y values, for regression models';
        'CM','Confusion matrix for discrimination and classification models (true classes in rows and predicted classes in column)';
        'PA','Prediction accuracy for binary discrimination models';
        'MCC','Matthews Correlation Coefficient for binary discrimination models';
        'layer', 'Whether this validation was an inner (single) validation or outer (double) validation';
        'scheme', 'Validation scheme to use for the validation';
        'n', 'Number of validation repeats';
        'k', 'The number of patitions the data was divided into';
        'permutation', 'Whether permutation was used or not';
        'confidence', 'Means, STD and confidence intervals for certain results';
        'optimization', 'Settings used for optimization of the model (inner validation'};
    Yp = validation;
end

end