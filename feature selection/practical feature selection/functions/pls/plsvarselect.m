function [model, selection, rank] = plsvarselect(model, method, varargin)

%DESCRIPTION:
%A function to optimize a PLS model trained with 'plstrain.m' by means of 
%predictor variable selection. Three methods for variable selection are 
%supported:
%- 1: Iteratively removing the predictor variable that has the lowest
%     absolute regression coefficient and re-calibrating the model on the
%     variables that are left. The selection of variables (with high
%     regression coefficient) leading to the model with the highest
%     performance is then considered optimal. This method is known as a
%     backpropagation selection based on regression coefficient. It is the
%     fastest method, but might produce suboptimal results.
%- 2: Also a backpropagation selection, but based on knock-out experiments. 
%     For each variable, it is determined how removing it would improve the 
%     model performance. The variable for which deletion improves the model 
%     most is then deleted and the knock-out experiment is repeated for the
%     variables that are left. This is repeated until no variables are
%     left. This method takes longer than method 1, especially when many
%     variables are present, but produces better results.
%- 3: A brute-forcing approach that tries all possible variable selections
%     using a full-factorial experimental design. This method is guarenteed
%     to find the optimal selection, but can take incredible long
%     especially when there are many predictor variables.
%
%This function is created only for PLS-R models and not PLS-DA models, and
%can only copy with single response variables. For multiple response
%variables, a response variable-specific selection of predictor variables
%is generally advised.
%
%INPUT:
%- model: The PLS model to perform the variable selection on, obtained 
%  using 'plstrain.m'.
%- method: The method used for optimization, see above.
%
%- The following name-value pairs can be used:
%  - 'optscheme': Validation scheme to use for the optimization of the
%    number of LVs:
%    - 'random' = Validation using a single random subset.
%    - 'kfold' = Random cross-validation using k folds.
%    - 'vblindsi' = Venetian blinds cross-validation on the sampling index
%      (recommended for time-series data).
%    - 'vblindsy'= Venetian blinds cross-validation on the first column of Y 
%      (default).
%    - 'leaveout' Leave-one-out cross-validation.
%    - 'cblock' = Contiguous block cross-validation.
%  - 'optn': Validation repeats to use (for non-deterministic schemes, 
%     default = 1).
%  - 'optk': For cross-validation, this is the number of partitions the data is
%    divided into (default = 5). For single test set validation, a fraction of 
%    1/k is taken out of the data for validation (so for k=5, 20% of the data 
%    is used for testing.
%  - 'optlv': Maximum number of LVs to consider for optimization of the
%    model.
%  - 'dontvalidate': Will make the algorithm optimize the variable
%    selection based on calibratd model performance rather than validated
%    performance. This will save a lot of time, but may cause an overfitted
%    selection.
%
%OUTPUT:
%- model: A copy of the model, but with the selected variables. This model
%  can be used for subsequent prediction, bootstrapping or validation. Only
%  when the estimated optimzation time is requested, this output will be
%  that estimated time.
%- selection: The selected variables (rows correspond to the indices of the 
%  original predictor data).
%- rank: For methods 1 and 2, these are the importance ranks for each of 
%  the variables. If the top number is 5, it means that the first variable
%  is the 5th most important. Method 3 does not rank on importance, and for 
%  this method rank is the experimental design matrix used.
%- performance: For methods 1 and 2, these are the RMSEs for the models
%  iteratively found while knocking out the variables. So, the RMSE at
%  position 10 of this vector correponds to the performance of the model
%  with the 10 highest-ranked variables. For method 3, this is the RMSE
%  found for each experiment.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 2021
%
%SYNTAX:
%[model, selection, rank, performance] = plsvarselect(model, 1, ...)
%[model, selection, rank, performance] = plsvarselect(model, 1, 'optscheme', 'vblinds', 'optn', 30, 'optk', '5')

%Parse input:
if nargin<2
    method=2;
end
if isfield(model, 'optimization')
    optscheme = model.optimization.optscheme;
    optn = model.optimization.optn;
    optk = model.optimization.optk;
    optlv = -model.optimization.LVmax;
else
    optscheme = 'vblinds';
    optn = 1;
    optk = 5;
    optlv = -model.LV;
end
dontvalidate = false;
for i=1:2:length(varargin)
    eval([varargin{i} ' = varargin{i+1};']);
end

%Optimize using method 1:
if method==1
	temp = model;
    performance = NaN(size(temp.X, 2), 1);
    rank = NaN(size(temp.X, 2), 1);
    indices = 1:size(temp.X, 2);
    while size(temp.X, 2)>0
        if dontvalidate
            performance(find(isnan(performance), 1, 'last')) = temp.calibration.RMSE;
        else
            temp.validation = plsvalidate(temp);
            performance(find(isnan(performance), 1, 'last')) = temp.validation.RMSE;
        end
        [~, i] = min(sum(abs(temp.B(2:end, 1)), 2));
        rank(find(isnan(rank), 1, 'last')) = indices(i);
        temp.X(:, i) = [];
        indices(i) = [];
        if size(temp.X, 2)>0
            temp = plstrain(temp.X, temp.Y, temp.Xpp, temp.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);
        end
    end
    [~, i] = min(performance);
    selection = rank(1:i);
    rank(rank) = 1:length(rank);
    model = plstrain(model.X(:, selection), model.Y, model.Xpp, model.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);

%Optimize using method 2:
elseif method==2
    temp = model;
    performance = NaN(size(temp.X, 2), 1);
    rank = NaN(size(temp.X, 2), 1);
    indices = 1:size(temp.X, 2);
    while size(temp.X, 2)>1
        if dontvalidate
            performance(find(isnan(performance), 1, 'last')) = temp.calibration.RMSE;
        else
            temp.validation = plsvalidate(temp);
            performance(find(isnan(performance), 1, 'last')) = temp.validation.RMSE;
        end
        P = NaN(size(temp.X, 2), 1);
        for j=1:length(P)
            M = plstrain(temp.X(:, setdiff(1:length(P), j)), temp.Y, temp.Xpp, temp.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);
            if dontvalidate
                P(j) = M.calibration.RMSE;
            else
                M.validation = plsvalidate(M);
                P(j) = M.validation.RMSE;
            end
        end
        [~, i] = min(P);
        rank(find(isnan(rank), 1, 'last')) = indices(i);
        temp.X(:, i) = [];
        indices(i) = [];
        if size(temp.X, 2)>0
            temp = plstrain(temp.X, temp.Y, temp.Xpp, temp.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);
        end
    end
    temp.validation = plsvalidate(temp);
    performance(find(isnan(performance), 1, 'last')) = temp.validation.RMSE;
    rank(1) = indices;
    [~, i] = min(performance);
    selection = rank(1:i);
    rank(rank) = 1:length(rank);
    model = plstrain(model.X(:, selection), model.Y, model.Xpp, model.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);

%Optimize using method 3:
elseif method==3
    doe = flipud(logical(fullfact(repmat(2, 1, size(model.X, 2)))-1));
    doe(end, :) = [];
    rmse = NaN(size(doe, 1), 1);
    for i=1:size(doe, 1)
        temp = plstrain(model.X(:, doe(i, :)), model.Y, model.Xpp, model.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);
        if dontvalidate
            rmse(i, 1) = temp.calibration.RMSE;
        else
            temp.validation = plsvalidate(temp);
            rmse(i, 1) = temp.validation.RMSE;
        end        
    end
    [~, i] = nanmin(rmse(:, 1));
    rank = doe;
    selection = find(rank(i, :));
    model = plstrain(model.X(:, selection), model.Y, model.Xpp, model.Ypp, optlv, 'optscheme', optscheme, 'optn', optn, 'optk', optk);
end

end