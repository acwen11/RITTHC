function lgltest(N)

[x,w,P,D,q,R] = lglgrid(N);

x = double(x);
w = double(w);
P = double(P);
D = double(D);
q = double(q);

xamr = ones(2*(N+1), 1);
xamr(1:N+1) = 0.5*x - 0.5;
xamr(N+2:2*N+2) = 0.5*x + 0.5;

% Check that the points are in the expected order
seq_err = max(abs(x-sort(x)))

fpoly = rand(N+1, 1);
fcoll = polyval(fpoly, x);

dfpoly = polyder(fpoly);
dfcoll = polyval(dfpoly, x);

fcollamr = polyval(fpoly, xamr);

% Relative maximum error in the derivative
diff_err = max(abs(D*fcoll-dfcoll))

ifpoly = polyint(fpoly')';
integral = polyval(ifpoly, 1) - polyval(ifpoly, -1);

% Quadrature error
quad_err = max(abs(w'*fcoll - integral))

% Prolongation error
prol_err = max(abs(R*fcoll - fcollamr))

% Restriction error
Rm = R \ eye(2*N+2);
restr_err = max(abs(Rm*R*fcoll - fcoll))
