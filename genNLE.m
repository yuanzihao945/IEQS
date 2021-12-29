function [nlef, nlefd, Xtrue] = genNLE(n)
% GENNLE - Generation of Nonlinear Functional Equations
% Input :
%   n       -  Dimensions * Functions (n * n)
%
% Output:  
%   nlef    -  n equations of nonlinear funcions
%   nlefd   -  diff of nonlinear funcions
%   Xtrue   -  true root of groups
%
% Usage:
%   [nlef, nlefd, Xtrue] = GENNLE(n, nset, nconset)
%
% See also INLEQS, ILEQS

% Author  : ZH.Yuan
% Update  : 2021/12/26 (First Version: 2021/12/26)
% Email   : zihaoyuan@whut.edu.cn (If any suggestions or questions)

Amethod = 3;
switch Amethod
    case 1
        A = n * eye(n) + rand(n);
    case 2
        A = 2 * eye(n) + rand(n);
        A = diag(diag(A, -1), -1) + diag(diag(A, 0), 0) + diag(diag(A, 1), 1);
    case 3
        A = n * eye(n) + rand(n);
        A = A' + A;
    case 4
        M = diag(rand(n, 1));
        Z = orth(rand(n, n));
        A = Z' * M * Z;
end

Xtrue = ones(n, 1) + rand(n, 1);
b = A * Xtrue;
k = randperm(14, 2);
[rNLEQ1, rNLEQD1] = TestFun(n, k(1));
[rNLEQ2, rNLEQD2] = TestFun(n, k(2));
nlef = @(x) [rNLEQ1(x) - rNLEQ1(Xtrue); rNLEQ2(x) - rNLEQ2(Xtrue); ...
    A(3 : end, :) * reshape(x, n, 1) - b(3 : end, :)];
nlefd = @(x) [rNLEQD1(x)'; rNLEQD2(x)'; A(3 : end, :)];
end


%% Generate Test Function
function [NLEQ, NLEQD] = TestFun(n, k)
% Select K n-ary nonlinear functions Randomly
if ~exist('k', 'var') || isempty(k)
    k = ceil(14 * rand);
end

A = [2; 2; zeros(n - 2, 1)] + rand(n, 1);
switch k
    case 1
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^2;
        NLEQD = @(x) 2 * sum(A .* reshape(x, n, 1)) * A;
    case 2
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^2);
        NLEQD = @(x) 2 * A .* reshape(x, n, 1);
    case 3
        NLEQ = @(x) log(sum(A .* reshape(x, n, 1)));
        NLEQD = @(x) A ./ sum(A .* reshape(x, n, 1));
    case 4
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^(3/2);
        NLEQD = @(x) (3/2) * sum(A .* reshape(x, n, 1))^(1/2) * A;
    case 5
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^(3/2));
        NLEQD = @(x) (3/2) * A .* reshape(x, n, 1).^(3/2);
    case 6
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^3;
        NLEQD = @(x) 3 * sum(A .* reshape(x, n, 1))^2 * A;
    case 7
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^(5/2);
        NLEQD = @(x) (5/2) * sum(A .* reshape(x, n, 1))^(3/2) * A;
    case 8
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^4;
        NLEQD = @(x) 4 * sum(A .* reshape(x, n, 1))^3 * A;
    case 9
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^3);
        NLEQD = @(x) 3 * A .* reshape(x, n, 1).^2;
    case 10
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^(5/2));
        NLEQD = @(x) (5/2) * A .* reshape(x, n, 1).^(3/2);
    case 11
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^3);
        NLEQD = @(x) 3 * A .* reshape(x, n, 1).^2;
    case 12
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^4);
        NLEQD = @(x) 4 * A .* reshape(x, n, 1).^3;
    case 13
        NLEQ = @(x) sum(A .* reshape(x, n, 1).^(7/2));
        NLEQD = @(x) (7/2) * A .* reshape(x, n, 1).^(5/2);
    case 14
        NLEQ = @(x) (sum(A .* reshape(x, n, 1)))^(7/2);
        NLEQD = @(x) (7/2) * sum(A .* reshape(x, n, 1))^(5/2) * A;
end

end