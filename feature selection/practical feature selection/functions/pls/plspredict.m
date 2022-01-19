function [Yp, RMSE, V] = plspredict(B, X, Yr, type)

%DESCRIPTION:
%A function to project existing data into an existing PLS(-DA) model. The
%function accepts model structures created with 'plstrain.m', or just 
%regression vectors. Be aware that in case of the latter, the function
%performs no preprocessing. A reference Y can be added to make the function
%calculate the prediction performance, which is useful for model
%validation. The function can return predicted values and prediction errors, 
%or a prediction-structure including all details of the projection.
%function can be used for predictions on models of the types 'regression',
%'discrimination' or 'classification'. The model type is auto-detected
%based on primarily on the model-structure and secondarily on the response
%references, but can be overridden with the input argument 'type'.
%
%INPUT:
%- B: The regression vector of the existing model.
%- X: The independent data to model.
%- Yr: The reference dependent data.
%- 'type': The model type to perform the prediction for.
%
%OUTPUT:
%- Yp: The predicted response values or classes.
%- RMSE: The performance of the model on the prediction data, regardless of
%  the model type, and only calcuted if the response reference is given. 
%  For discrimination and classification, this number is calculated before 
%  rounding the responses to a dummy-matrix.
%- V: A 2-by-1 vector giving the fraction of variance explained by the entire
%  model. The first number is the variance in X, the second one is the
%  variance in Y. For multiple Y variables, this represents the variance 
%  of Y in the entire Y-block. The variance in X is only calculated when a 
%  model-structure is inputted; the variance in Y is only calculated when a
%  Yr is inputted.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 2021
%
%SYNTAX:
%[Yp, RMSE, V] = plspredict(B, X, Yr, type)
%[prediction] = plspredict(model, X, [], 'regression')

%Handle missing input:
if nargin<4
    type = [];
end
if nargin<3
    Yr = [];
end

%Intitialize prediction-structure:
Yp = [];
RMSE = [];
V = [];
prediction.type = [];
prediction.Yp = [];
prediction.RMSE = [];
prediction.V = [];
prediction.X = X;
if exist('Yr')==1
    prediction.Yr = Yr;
end
    
%Auto-select the correct PLS model type (if none is given):
if isempty(type)
    if isstruct(B)
        type = B.type;
    elseif ~isempty(Yr)
        if isstring(Yr)
            type = 'discrimination';
        elseif sum(round(Yr(:))-Yr(:))==0
            if size(Yr, 2)==1
                type = 'discrimination';
            elseif ~any(sum(Yr, 2)~=1)
                type = 'discrimination';
            else
                type = 'classification';
            end
        end
    else
        type = 'regression';
    end
end

%Perform prediction for regression vectors:
if isnumeric(B)
    Yp = (X * B(2:end, :)) + B(1, :);
    T2p = [];
    Qp = [];
    
%Perform prediction for model structures, uncluding preprocessing:
elseif isstruct(B)
    
    %Preprocess X:
    if B.Xpp>0
        X = X - (ones(size(X, 1), 1) * mean(B.X));
    end
    if B.Xpp>1
        X = X ./ (ones(size(X, 1), 1) * std(B.X));
    end
    
    %Replace Infs and NaNs by zeros:
    X(isinf(X)) = 0;
    X(isnan(X)) = 0;
    
    %Predict:
    Yp = (X * B.B(2:end, :)) + B.B(1, :);
    Tp = (X * B.P);
    T2p = sum((Tp ./ (ones(size(Tp, 1), 1) * (B.V(1, :).^0.5))).^2, 2);
    Qp = sum((X - Tp*B.P').^2, 2);
   
    %Preprocess Y:
	if B.Ypp>1
        Yp = Yp .* (ones(size(Yp, 1), 1) * std(B.Y));
	end
    if B.Ypp>0
    	Yp = Yp + (ones(size(Yp, 1), 1) * mean(B.Y));
    end
    
    %Calculate explained variance in X:
    Xf = B.T * B.P';
    V(1, 1) = sum(Xf(:).^2) / sum(X(:).^2);
    
end

%Transform Yr to dummy matrix for discriminations and classifications:
if (~isempty(Yr)) && (strcmp(type, 'discrimination') || strcmp(type, 'classification'))
    if size(Yr, 2)==1
        if isstruct(B)
            if isfield(B, 'Yclasses')
                try
                    Yr = double(repmat(B.Yclasses', size(Yr, 1), 1) == Yr);
                catch
                    warning('Format of reference responses does not match format of response in model-structure. Ignoring responses in model-structure. This can causes problems if not all classes in the model-structure are in the reference responses.');
                    [Yclasses, ~, Yr] = unique(Yr, 'stable');
                    Yr = double(repmat(unique(Yr)', size(Yr, 1), 1) == repmat(Yr, 1, length(unique(Yr))));
                end
            else
                [Yclasses, ~, Yr] = unique(Yr, 'stable');
                Yr = double(repmat(unique(Yr)', size(Yr, 1), 1) == repmat(Yr, 1, length(unique(Yr))));
            end
        else
            [Yclasses, ~, Yr] = unique(Yr, 'stable');
            Yr = double(repmat(unique(Yr)', size(Yr, 1), 1) == repmat(Yr, 1, length(unique(Yr))));
        end
    else
        Yr = double(Yr-min(Yr(:))>0.5);
    end
    prediction.Yr = Yr;
end

%Calculate model RMSE and explained variance in Y:
if (~isempty(Yr))
    RMSE = sqrt(mean((Yr(:)-Yp(:)).^2));
    if isstruct(B)
        if B.Ypp>0
            Yr = Yr - (ones(size(Yr, 1), 1) * mean(B.Y));
        end
        if B.Ypp>1
            Yr = Yr ./ (ones(size(Yr, 1), 1) * std(B.Y));
        end
        Yf = (X * B.B(2:end, :)) + B.B(1, :);
    else
        Yf = (X * B(2:end, :)) + B(1, :);
    end
    V(2, 1) = sum(Yf(:).^2) / sum(Yr(:).^2);
end

%Calculate R for regression models:
if (~isempty(Yr)) && strcmp(type, 'regression')
    R = min(min(corrcoef(Yr(:), Yp(:))));
end

%Transform Yp to dummy matrix in case of discrimination or classification:
if strcmp(type, 'discrimination')
    Yp = double(Yp == repmat(max(Yp, [], 2), 1, size(Yp, 2)));
elseif strcmp(type, 'classification')
    Yp = double(Yp>=0.5);
end

%Calculate CM for discrimination and classification models (true classes in rows and predicted classes in column), and PA and MCC if possible:
if (~isempty(Yr)) && (strcmp(type, 'discrimination') || strcmp(type, 'classification'))
    CM = zeros(size(Yr, 2));
    for i=1:size(CM, 1)
        CM(i, :) = sum(repmat(Yr(:, i), 1, size(Yr, 2)) .* Yp);
    end
    if (size(Yr, 2)==1) && (~any(sum(Yr, 2)~=1))
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

%Prepare prediction-structure, if requested to do so:
if nargout==1
    prediction.Yp = Yp;
    prediction.RMSE = RMSE;
    prediction.V = V;
    prediction.type = type;
    if exist('R')==1
        prediction.R = R;
    end
    if exist('CM')==1
        prediction.CM = CM;
    end
    if exist('PA')==1
        prediction.PA = PA;
    end
    if exist('MCC')==1
        prediction.MCC = MCC;
    end
    if ~isempty('Tp')
        prediction.Tp = Tp;
    end
    if ~isempty('T2p')
        prediction.Tp = Tp;
    end
    if ~isempty('Qp')
        prediction.Tp = Tp;
    end
    prediction.info = {
        'type','Prediction type (regression, discrimination or classification)';
        'Yp','Predicted response values or classes';
        'X','Inputted independent data';
        'Yr','Reference dependent data';
        'RMSE','Root mean squared error of prediction';
        'V','Fraction of explained variance in X and Y by each latent variable';
        'R','Correlation between prediction and reference Y values, for regression models';
        'CM','Confusion matrix for discrimination and classification models (true classes in rows and predicted classes in column)';
        'PA','Prediction accuracy for binary discrimination models';
        'MCC','Matthews Correlation Coefficient for binary discrimination models';
        'Tp','Projected scores of X in model';
        'T2p','Hotellings T-squared values for each sample';
        'Qp','Q-residual for each sample'};
    Yp = prediction;
end

end