% Computes the Gauss-Legendre-Lobatto derivative matrix using a formula from
% Canuto, Hussaini, Quarteroni, Zang - Spectral Methods (2006)
%
% x : abscissas (assumed to be sorted)
% w : weights
% P : Vandermonde matrix of the Legendre polynomials
%
% D : discrete derivative operator: u' ~ D*u
function D = lgldiff(x, w, P)

N1 = length(x);
N  = N1 - 1;

D = mp(zeros(N1, N1));

D(1,1)   = - (N + 1.0)*N /4.0;
D(N1,N1) = - D(1,1);

for i = 1:N1
    for j = 1:N1
        if i ~= j
            D(i,j) = P(i,N1)/P(j,N1) * 1.0/(x(i)-x(j));
        end
    end
end
