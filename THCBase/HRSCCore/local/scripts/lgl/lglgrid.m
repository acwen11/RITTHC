% Computes the collocation points and weights of the Guass-Legendre-Lobatto
% quadrature on [-1, 1].
%
% x  : abscissas
% w  : weights
% P  : Vandermonde matrix of the Legendre polynomials
% D  : discrete derivative operator
% q  : interpolation weights: q_k = 1 / \Pi_{i\neq k} (x_k - x_i)
% R  : prolongation operator for fixed factor 2 AMR
%
% Rm : pseudo-inverse of R
%
% N : number of grid points
function [x, w, P, D, q, R] = lglgrid(N)

[x, w, P] = lglnodes(N);
D = lgldiff(x, w, P);

N1 = N+1;
q = mp(ones(N1,1));
for i = 1:N1
    for j = 1:N1
        if i ~= j
            q(i) = q(i) * (x(i) - x(j));
        end
    end
end
q = 1.0 ./ q;

% We split the element [-1,1] into two identical non-overlapping elements
% [-1,0] and [0,1]
R = ones(2*N1,N1);

xamr = ones(1, 2*N1);
xamr(1:N1) = 0.5*x - 0.5;
xamr(N1+1:2*N1) = 0.5*x + 0.5;

% Compute Lagrange basis
l = (vander(double(x))\eye(N1))';
for i = 1:N1
    R(:,i) = polyval(l(i,:), xamr);
end
