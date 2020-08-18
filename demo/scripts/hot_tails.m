clear all

%------------------------------------------------------------------------------|
%   Set plasma parameters
%------------------------------------------------------------------------------|
time    = 0:5e-6:2e-3;
t_star  = 1.5e-4*ones(1,1);
n0      = 3.0e19*ones(1,1);
nmax    = 15.0e19*ones(1,1);
T0      = 7.0e3*ones(1,1);
Tf      = 10.0*ones(1,1);
Epar    = (0.01 - 1.)*exp(-time/5.e-4) + 1.;
T       = Tf(1,1) + (T0(1,1)-Tf(1,1)).*exp(-(time/t_star));
n       = nmax(1,1) + (n0(1,1)-nmax(1,1)).*exp(-(time/t_star));

%------------------------------------------------------------------------------|
%   Initialize hot-tails
%------------------------------------------------------------------------------|
H = hot_tails_class('t_star',t_star,'n0',n0,'nmax',nmax,'T0',T0,'Tf',Tf);

%------------------------------------------------------------------------------|
%   Calculate evolution of hot-tail population
%------------------------------------------------------------------------------|
for i = 1:length(time)
    Ht(i,:) = H.calculate(time(i),Epar(i)*ones(1,1));
end

%------------------------------------------------------------------------------|
%   Write results to file
%------------------------------------------------------------------------------|
fileID = fopen('../dat/hot_tails_matlab.dat', 'w');
fprintf(fileID,'%-20s', '# Time (s)', '  n_hot (m**-3)', '  n_e (m**-3)', ...
    '  T_e (ev)', '  E_par (V/m)', '  v_c (v_th0)', '  tau');
fprintf(fileID,'\n');
for i = 1:length(time)
    fprintf(fileID,'%19.12e ', time(i), Ht(i), n(i), T(i), Epar(i), 0., 0.);
    fprintf(fileID,'\n');
end
fclose(fileID);

%------------------------------------------------------------------------------|