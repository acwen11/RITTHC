% Generates a sequence of LGL Grids and outputs them to file

% Minimum number of points
istart = 1;
% Maximum number of points
iend   = 25;
% Precision
iprec  = 20;

for N = istart:iend
    [x, w, P, D, q, R] = lglgrid(N);

    Q = inv(P);

    x = double(x);
    w = double(w);
    P = double(P);
    Q = double(Q);
    D = double(D);
    q = double(q);

    Rm = R \ eye(2*N+2);

    fname = strcat('nodes.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    for i = 1:length(x)-1
        fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), x(i));
    end
    fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), x(length(x)));
    fclose(fid);

    fname = strcat('weights.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    for i = 1:length(w)-1
        fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), w(i));
    end
    fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), w(length(w)));
    fclose(fid);

    fname = strcat('icoeff.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    for i = 1:length(q)-1
        fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), q(i));
    end
    fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), q(length(q)));
    fclose(fid);

    fname = strcat('diffop.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    sz = size(D);
    for i = 1:sz(1)
        for j = 1:sz(2)
            if i == sz(1) & j == sz(2)
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), D(i,j));
            else
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), D(i,j));
            end
        end
    end
    fclose(fid);

    fname = strcat('dltop.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    sz = size(Q);
    for i = 1:sz(1)
        for j = 1:sz(2)
            if i == sz(1) & j == sz(2)
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), Q(i,j));
            else
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), Q(i,j));
            end
        end
    end
    fclose(fid);

    fname = strcat('idltop.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    sz = size(P);
    for i = 1:sz(1)
        for j = 1:sz(2)
            if i == sz(1) & j == sz(2)
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), P(i,j));
            else
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), P(i,j));
            end
        end
    end
    fclose(fid);

    fname = strcat('prolongation.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    sz = size(R);
    for i = 1:sz(1)
        for j = 1:sz(2)
            if i == sz(1) & j == sz(2)
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), R(i,j));
            else
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), R(i,j));
            end
        end
    end
    fclose(fid);

    fname = strcat('restriction.', int2str(N), '.inc');
    fid = fopen(fname, 'w');
    sz = size(Rm);
    for i = 1:sz(1)
        for j = 1:sz(2)
            if i == sz(1) & j == sz(2)
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g\n'), Rm(i,j));
            else
                fprintf(fid, strcat('%1.', mat2str(iprec), 'g,\n'), Rm(i,j));
            end
        end
    end
    fclose(fid);
end
