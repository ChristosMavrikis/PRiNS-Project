function [B, T, U, P, Q, V] = plstrain(X, Y, Xpp, Ypp, LV, varargin)

%DESCRIPTION:
%A function to perform Partial Least Squares modelling. You can specify the
%number of latent variables to model, or let the function optimize this
%automatically using a tweakable validation scheme. You can ask the
%individual output arguments of the function, or you can ask to return
%everything in a model-structure. This model-structure contains all output
%arguments given below, plus some additional information.
%
%The models can be of the type 'regression', 'discrimination' and
%'classification'. For discrimination, samples belong to either one of two
%classes. For classification, samples may belong to multiple or to no
%classes. The function auto-selects the correct model based on the input,
%but you can also force the model to a certain type (see below). For
%discrimination and classification, the input is always transformed to a
%dummy-matrix.
%
%INPUT:
%- X: The independent data to model.
%- Y: The dependent data to model.
%- Xpp: The preprocessing method for the independent data (0=none, 1=mean
%  centering, 2=autoscaling).
%- Ypp: The preprocessing method for the dependent data (0=none, 1=mean
%  centering, 2=autoscaling).
%- LV: The number of latent variables to model:
%  - LV = n: Model n latent variables.
%  - LV = Inf: Model maximum number of latent variables.
%  - LV = -n: Optimize number of latent variables between 1 and n using
%    cross-validation.
%
%- The following name-value pairs can be used:
%  - 'type': Type of PLS-model ('regression', 'discrimination' or
%    'classification'.
%  - 'optscheme': Validation scheme to use for the optimization of the
%    number of LVs:
%    - 'random': Validation using a single random subset.
%    - 'kfold': Random cross-validation using k folds.
%    - 'vblindsi': Venetian blinds cross-validation on the sampling index
%      (recommended for time-series data).
%    - 'vblindsy': Venetian blinds cross-validation on the first column of Y 
%      (default).
%    - 'leaveout': Leave-one-out cross-validation.
%    - 'cblock': Contiguous block cross-validation.
%  - 'optn': Validation repeats to use (for non-deterministic schemes, 
%     default = 1).
%  - 'optk': For cross-validation, this is the number of partitions the data is
%    divided into (default = 5). For single test set validation, a fraction of 
%    1/k is taken out of the data for validation (so for k=5, 20% of the data 
%    is used for testing.
%  - 'signstable': Whether to use the signstable SVD by Rasmus Bro 
%    (https://doi.org/10.1002/cem.1122) or not (default).
%  - 'algorithm': The algorithm to use for model calculation: 'NIPALS'
%    (default) or 'SIMPLS'.
%
%OUTPUT:
%- B: Regression vector. The first element is B0.
%- T: X-scores.
%- U: Y-scores.
%- P: X-loadings.
%- Q: Y-loadings.
%- V: Fraction of variance explained by each Latent Variable. The columns
%  are the Latent Variables, the first row the variation in X and the second
%  row the variation in Y. For multiple Y variables, this represents the
%  variance of Y in the entire Y-block.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 2021
%
%SYNTAX:
%[B, T, U, P, Q, V] = plstrain(X, Y, Xpp, Ypp, LV, 'optscheme', 'kfold', 'n', 10, 'k', 5)
%[model] = plstrain(X, Y, Xpp, Ypp, LV, ...)

%Parse input:
type = [];
optscheme = 'vblinds';
optn = 1;
optk = 5;
signstable = false;
algorithm = 'NIPALS';
for i=1:2:length(varargin)
    eval([varargin{i} ' = varargin{i+1};']);
end

%Correct input:
if optk<1
    optk = ceil(1/optk);
end

%Initialize the model-structure and save the raw X and Y data. The original raw 
%data is used by the algorithm later on.
model.B = [];
model.T = [];
model.U = [];
model.P = [];
model.Q = [];
model.V =[];
model.LV = [];
model.X = X;
model.Y = [];

%Check and correct the number of latent variables requested for modelling:
if abs(LV)>min(size(X))
    LV = sign(LV) * min(size(X));
end
if size(X, 2) == 1
    LV = 1;
end

%Auto-select the correct PLS model type (if none is given):
if isempty(type)
    if isstring(Y)
        type = 'discrimination';
    elseif sum(round(Y(:))-Y(:))==0
        if size(Y, 2)==1
            type = 'discrimination';
        elseif ~any(sum(Y, 2)~=1)
            type = 'discrimination';
        else
            type = 'classification';
        end
    else
        type = 'regression';
    end
end

%Transform the response in case of discrimination or classification:
Yclasses = [];
if strcmp(type, 'discrimination')
    if size(Y, 2)==1
        [Yclasses, ~, Y] = unique(Y, 'stable');
        Y = double(repmat(unique(Y)', size(Y, 1), 1) == repmat(Y, 1, length(unique(Y))));
    else
        Y = Y - min(Y(:));
    end
    Ypp = 0;
elseif strcmp(type, 'classification')
    Y = Y - min(Y(:));
    Ypp = 0;
end
model.Y = Y;

%Optimize the number of latent variables using (cross-)validation. The
%settings to use for cross-validation can be set by the user, but on
%default a Venetian blinds on the sorted Y-variable is used.
if LV<0
    
    %Copy maximum number of LVs to model-structure:
    model.optimization.LVmax = abs(LV);
    model.optimization.optscheme = optscheme;
    model.optimization.optn = optn;
    model.optimization.optk = optk;
    model.optimization.Ypcal = [];
    model.optimization.Ypval = [];
    
    %Determine validation set(s):
    if strcmp(optscheme, 'vblindsi') || strcmp(optscheme, 'vblinds')
        optn = 1;
        test = false(optn, optk, size(model.Y, 1));
        for i=1:size(test, 2)
            test(1, i, i:size(test, 2):size(test, 3)) = true;
        end
    elseif strcmp(optscheme, 'vblindsy')
        optn = 1;
        test = false(optn, optk, size(model.Y, 1));
        [~, s] = sort(model.Y(:, 1));
        for i=1:size(test, 2)
            test(1, i, s(i:size(test, 2):size(test, 3))) = true;
        end
    elseif strcmp(optscheme, 'kfold')
        test = false(optn, optk, size(model.Y, 1));
        for n = 1:size(test, 1)
            s = randperm(size(test, 3));
            for i=1:size(test, 2)
                test(n, i, s(i:size(test, 2):size(test, 3))) = true;
            end
        end
    elseif strcmp(optscheme, 'random')
        test = false(optn, 1, size(model.Y, 1));
        for n = 1:size(test, 1)
            s = randperm(size(test, 3));
            test(n, 1, s(1:floor(size(test, 3) * (1/optk)))) = true;
        end
    elseif strcmp(optscheme, 'leaveout')
        optn = 1;
        test = false(optn, optk, size(model.Y, 1));
        optk = size(model.Y, 1);
        test = false(1, size(model.Y, 1), size(model.Y, 1));
        test(1, find(eye(size(model.Y, 1)))) = true;
    elseif strcmp(optscheme, 'cblock')
        optn = 1;
        test = false(optn, optk, size(model.Y, 1));
        for i=1:size(test, 2)
            test(optn, i, floor((i-1) * (size(test, 3) / size(test, 2)) + 1):floor((i) * (size(test, 3) / size(test, 2)))) = true;
        end
	else
        error('Invalid validation scheme. Please check the help-information and your input.');
    end
    
    %Initialize validation output:
    errors = NaN([size(test, 1), abs(LV), size(Y)]);
    Ypcal = NaN([size(test, 1), size(test, 2), abs(LV), size(Y)]);
    Ypval = NaN([size(test, 1), abs(LV), size(Y)]);

    %For each repeat:
    for n = 1:size(test, 1)
    
        %For each fold:
        for k=1:size(test, 2)

            %Divide test and training set:
            Xtrain = X(~test(n, k, :), :);
            Ytrain = Y(~test(n, k, :), :);
            Xtest = X(test(n, k, :), :);
            Ytest = Y(test(n, k, :), :);

            %Preprocess X of train and test set:
            if Xpp>0
                Xtest = Xtest - (ones(size(Xtest, 1), 1) * mean(Xtrain));
                Xtrain = Xtrain - (ones(size(Xtrain, 1), 1) * mean(Xtrain));
            end
            if Xpp>1
                Xtest = Xtest ./ (ones(size(Xtest, 1), 1) * std(Xtrain));
                Xtrain = Xtrain ./ (ones(size(Xtrain, 1), 1) * std(Xtrain));
            end

            %For each latent variable:
            for ilv = 1:abs(LV)

                %Calibrate, test and save performance:
                B = pls(Xtrain, Ytrain, Xpp, Ypp, ilv, signstable, algorithm);
                
                %Calibration performance:
                Yp = (Xtrain * B(2:end, :)) + B(1, :);
                if Ypp>1
                    Yp = Yp .* (ones(size(Yp, 1), 1) * std(Ytrain));
                end
                if Ypp>0
                    Yp = Yp + (ones(size(Yp, 1), 1) * mean(Ytrain));
                end
                Ypcal(n, k, ilv, ~test(n, k, :), :) = Yp;
                
                %Cross-validation performance:
                Yp = (Xtest * B(2:end, :)) + B(1, :);
                if Ypp>1
                    Yp = Yp .* (ones(size(Yp, 1), 1) * std(Ytrain));
                end
                if Ypp>0
                    Yp = Yp + (ones(size(Yp, 1), 1) * mean(Ytrain));
                end
                Ypval(n, ilv, test(n, k, :), :) = Yp;
                errors(n, ilv, test(n, k, :), :) = Ytest - Yp;
            end
        end
    end
    
	%Select the optimal number of components:
    errors = squeeze(nanmean(errors.^2, 1).^0.5);
    errors = squeeze(nanmean(errors(:, :).^2, 2).^0.5);
    [~, LV] = min(errors);
    
    %Average results:
    Ypcal = nanmean(Ypcal, 1);
    Ypcal = squeeze(nanmean(Ypcal, 2));
    Ypcal = squeeze(Ypcal(LV, :, :));
    Ypcal = reshape(Ypcal, size(Y));
    model.optimization.Ypcal = Ypcal;
    Ypval = squeeze(nanmean(Ypval, 1));
    Ypval = squeeze(Ypval(LV, :, :));
    Ypval = reshape(Ypval, size(Y));
    model.optimization.Ypval = Ypval;
    
end

%Train (optimized) model:
[B, T, U, P, Q, V, signstable, algorithm, W] = pls(X, Y, Xpp, Ypp, LV, signstable, algorithm);

%Finalize the model-type output, if requested:
if nargout==1
    model.B = B;
    model.T = T;
    model.U = U;
    model.P = P;
    model.Q = Q;
    model.V = V;
    model.W = W;
    model.LV = LV;
    if ~isempty(Yclasses)
        model.Yclasses = Yclasses;
    end
    model.Xpp = Xpp;
    model.Ypp = Ypp;
    model.type = type;
    model.signstable = signstable;
    model.algorithm = algorithm;
    model.info = {
        'B','Regression vector';
        'T','X scores of each sample';
        'U','Y scores of each sample';
        'P','X loadings';
        'Q','Y loadings';
        'V','Fraction of explained variance in X and Y for each latent variable';
        'LV','Number of latent variables fitted';
        'X','Input independent data';
        'Y','Input dependent data';
        'Yclasses','The original names for the classes (in case the model type is discrimination or classification)';
        'Xpp','Preprocessing method for X';
        'Ypp','Preprocessing method for Y';
        'type', 'Type of PLS-model (''regression'', ''discrimination'' or ''classification''';
        'signstable', 'Whether the signstable SVD was used or not.';
        'algorithm', 'Algorithm used to calculate the PLS model.';
        'calibration','Prediction of the training data in the model'};
    if isfield(model, 'optimization')
        model.info = [model.info; {
            'optmization.LVmax', 'Maximum number of LVs condisidered during optimization';
            'optmization.optscheme', 'Validation scheme used for the optimization of the number of LVs';
            'optmization.optn', 'Validation repeats used for the optimization of the number of LVs';
            'optmization.optk', 'For cross-validation, this is the number of partitions the data is divided into (default = 5). For single test set validation, a fraction of 1/k is taken out of the data for validation (so for k=5, 20% of the data is used as test set)'
            'optmization.Ypcal', 'Predicted Y-values for model calibration, can be used to calculate the calibration error.'
            'optmization.Ypval', 'Predicted Y-values for model optimization (validation), can be used to calculate the cross-validation error.'}];
    end
        
    %Predict the calibration data in the model:
    if exist('plspredict.m')
        model.calibration = plspredict(model, model.X, model.Y);
    end
    
    %Output model-structure as primary output:
    B = model;
    
end

end