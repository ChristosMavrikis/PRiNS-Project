function [B, T, U, P, Q, V, signstable, algorithm, W] = pls(X, Y, Xpp, Ypp, LV, signstable, algorithm)

%DESCRIPTION:
%A basic function to perform pls-modelling. You can use it to train models,
%but you can also use 'plstrain.m', which give you many more options,
%including an optimization routine. To calculate the model, you can choose
%between the NIPALS (default) or SIMPLS algorithms.
%
%INPUT:
%- X: The independent data to model.
%- Y: The dependent data to model.
%- Xpp: The preprocessing method for the independent data (0=none, 1=mean
%  centering, 2=autoscaling).
%- Ypp: The preprocessing method for the dependent data (0=none, 1=mean
%  centering, 2=autoscaling).
%- LV: The number of latent variables to model.
%- signstable: Whether to use the signstable SVD by Rasmus Bro 
%  (https://doi.org/10.1002/cem.1122) or not (default).
%- algorithm: The algorithm to use for model calculation: 'NIPALS'
%  (default) or 'SIMPLS'.
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
%- signstable: Whether the signstable SVD was used or not.
%- algorithm: The algorithm used to calculate this particular model.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 202!
%
%SYNTAX:
%[B, T, U, P, Q, V] = pls(X, Y, Xpp, Ypp, LV)
%[B, T, U, P, Q, V, signstable, algorithm] = pls(X, Y, Xpp, Ypp, LV, signstable, algorithm)

%Check input:
if nargin<7
    algorithm = 'NIPALS';
end
if nargin<6
    signstable = false;
end
if nargin<5
    LV = min(size(X));
end
if nargin<4
    Ypp = 0;
end
if nargin<3
    Xpp = 0;
end
if LV > min(size(X))
    LV = min(size(X));
end

%Check if the function for signstable SVD is present:
if signstable
    if ~exist('svd_signstable.m')
        signstable = false;
        warning('Function for signstable SVD could not be found (svd_signstable.m). Using regular SVD instead.');
    end    
end

%Turn of warnings:
warning('off');

%Preprocess data:
if Xpp>0
    X = X - (ones(size(X, 1), 1) * mean(X));
end
if Xpp>1
    X = X ./ (ones(size(X, 1), 1) * std(X));
end
if Ypp>0
    Y = Y - (ones(size(Y, 1), 1) * mean(Y));
end
if Ypp>1
    Y = Y ./ (ones(size(Y, 1), 1) * std(Y));
end

%Check for NaNs:
X(isinf(X)) = NaN;
Xtemp = X;
X(:, sum(isnan(X))>0) = [];

%Initialize outputs:
T = zeros(size(X, 1), LV);
U = zeros(size(X, 1), LV);
P = zeros(size(X, 2), LV);
Q = zeros(size(Y, 2), LV);
W = zeros(size(X, 2), LV);

%Calculate model using SIMPLS, if requested:
if strcmpi(algorithm, 'SIMPLS')
    algorithm = 'SIMPLS';
    
    %Calculate latent variables:
    V = zeros(size(X, 2), LV);
    S = X' * Y;
    for i=1:LV

        %Calculate orthonormal basis:
        if signstable
            [ri, si, ci] = svd_signstable(S, 'econ');
        else
            [ri, si, ci] = svd(S, 'econ');
        end
        ri = ri(:, 1);
        ci = ci(:, 1);
        si = si(1);

        %Calculate loadings:
        ti = X*ri;
        normti = norm(ti);
        ti = ti ./ normti;
        P(:, i) = X'*ti;
        qi = si*ci/normti;
        Q(:, i) = qi;

        %Calculate scores and weights:
        T(:, i) = ti;
        U(:, i) = Y*qi;
        W(:, i) = ri ./ normti;

        %Update orthonormal basis:
        vi = P(:, i);
        for repeat = 1:2
            for j=1:i-1
                vj = V(:, j);
                vi = vi - (vj'*vi)*vj;
            end
        end
        vi = vi ./ norm(vi);
        V(:, i) = vi;

        %Deflate covariance matrix:
        S = S - vi*(vi'*S);
        Vi = V(:, 1:i);
        S = S - Vi*(Vi'*S);
    end

    %Orthogonalize Y scores to preceeding X scores:
    for i=1:LV
        ui = U(:, i);
        for repeat = 1:2
            for j=1:i-1
                tj = T(:, j);
                ui = ui - (tj'*ui)*tj;
            end
        end
        U(:, i) = ui;
    end

    B = NaN;
    V = NaN;

    %Calculate regression vector:
    B = W * Q';
    B = [mean(Y, 1) - mean(X, 1)*B; B];

    %Calculate fraction of explained variance:
    V = [sum(abs(P).^2,1) ./ sum(sum(abs(X).^2,1)); sum(abs(Q).^2,1) ./ sum(sum(abs(Y).^2,1))];

%Calculate model using NIPALS, in all other cases:
else
    algorithm = 'NIPALS';
    
    %Build outer model:
    E = X;
    F = Y;
    for i = 1:LV
        S = E' * F;
        if signstable
            [UU, D, V] = svd_signstable(S, 'econ');
        else
            [UU, D, V] = svd(S, 'econ');
        end
        W(:, i) = UU(:, 1);
        T(:, i) = E * W(:, i);
        P(:, i) = E' * T(:, i) / (T(:, i)' * T(:, i));
        Q(:, i) = F' * T(:, i) / (T(:, i)' * T(:, i));
        U(:, i) = F * Q(:, i) / (Q(:, i)' * Q(:, i));
        E = E - T(:, i) * P(:, i)';
        F = F - T(:, i) * Q(:, i)';
    end

    %Calculate fractions of explained variance:
    V = zeros(2, LV);
    for i=1:LV
        Xf = T(:, i) * P(:, i)';
        V(1, i) = sum(Xf(:).^2) / sum(X(:).^2);
        D = W(:, 1:i) * inv(P(:, 1:i)' * W(:, 1:i)) * Q(:, 1:i)';
        B = [mean(mean(Y, 1) - (X * D)); D];
        Yf = X * B(2:end, :) + B(1, :);
        V(2, i) = sum(Yf(:).^2) / sum(Y(:).^2);
    end
    V(2, 2:end) = V(2, 2:end) - V(2, 1:end-1);

    %Build inner model:
    D = W * inv(P' * W) * Q';
    B = [mean(mean(Y, 1) - (X * D)); D];
end

%Correct for NaNs:
if size(Xtemp, 2) > size(X, 2)
    Btemp = zeros(size(Xtemp, 2)+1, size(B, 2));
    Ptemp = zeros(size(Xtemp, 2), size(P, 2));
    Btemp(find(sum(isnan(Xtemp))==0)+1, :) = B(2:end, :);
    Ptemp(find(sum(isnan(Xtemp))==0), :) = P;
    B = Btemp;
    P = Ptemp;
end

%Turn on warnings:
warning('on');
    
end