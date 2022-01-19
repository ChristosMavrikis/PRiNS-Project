function [U, D, V] = svd_signstable(X, economy)

%DESCRIPTION:
%This function assigns meaningful, deterministic signs to the output of a
%singular value decomposition, and provides a sign-stable manner for
%performing such decomposition. Based on Bro, Rasmus, Evrim Acar, and 
%Tamara G. Kolda. "Resolving the sign ambiguity in the singular value 
%decomposition." Journal of Chemometrics 22.2 (2008): 135-140.
%
%INPUT:
%- X: Data matrix to decompose.
%- econmy: Whether to quicker-to-calculate economy SVD (default) or not.
%
%OUTPUT:
%- U: Sign-corrected right-singular vectors of X.
%- D: Singular values of X.
%- V: Sign-corrected left-singular vectors of X.
%
%AUTHOR:
%Tim Offermans, Radboud University Nijmegen (The Netherlands), November 2021
%
%SYNTAX:
%[U, D, V] = svd_signstable(X, economy)

%Check input:
if nargin<2
    economy = true;
end
if strcmp(economy, 'econ')
    economy = true;
end

%Perform the (sign-unstable) SVD:
if economy
    [U, D, V] = svd(X, 'econ');
else
    [U, D, V] = svd(X);
end

%Step 1:
K = length(D);
s_left = NaN(1, K);
for k=1:K
    Y = X - U(:, setdiff(1:K, k)) * D(setdiff(1:K, k), setdiff(1:K, k)) * (V(:, setdiff(1:K, k))');
    s_left_parts = NaN(1, size(Y, 2));
    for j=1:size(Y, 2)
        temp_prod = (U(:, k)') * Y(:, j);
        s_left_parts(j) = sign(temp_prod) * temp_prod^2;
    end
    s_left(k) = sum(s_left_parts);
end

%Step 2:
s_right = NaN(1, K);
for k=1:K
    Y = X - U(:, setdiff(1:K, k)) * D(setdiff(1:K, k), setdiff(1:K, k)) * (V(:, setdiff(1:K, k))');
    s_right_parts = NaN(1, size(Y, 1));
    for i=1:size(Y, 1)
        temp_prod = (V(:, k)') * (Y(i, :)');
        s_right_parts(i) = sign(temp_prod) * temp_prod^2;
    end
    s_right(k) = sum(s_right_parts);
end

%Step 3:
for k=1:K
    if (s_right(k) * s_left(k)) < 0
        if s_left(k) < s_right(k)
            s_left(k) = -s_left(k);
        else
            s_right(k) = -s_right(k);
        end
    end
end
U = U * diag(sign(s_left));
V = V * diag(sign(s_right));

end