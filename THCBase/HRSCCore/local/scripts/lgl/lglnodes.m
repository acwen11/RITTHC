function [x,w,P]=lglnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Modified by David Radice - 07/14/2011 to reorder the nodes
%                          - 07/26/2011 to use the MP toolbox
% Contact: gregvw@chtm.unm.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=-mp(cos(pi*(0:N)/N)');

% The Legendre Vandermonde Matrix
P=mp(zeros(N1,N1));

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>2^(-100)

    xold=x;

    P(:,1)=1;    P(:,2)=x;

    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end

    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );

end

P(:,1)=1;
P(:,2)=x;
for k=2:N
    P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
end

w=2./(N*N1*P(:,N1).^2);
