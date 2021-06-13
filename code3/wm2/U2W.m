function [ W ] = U2W( U )
% transform U to W

global gamma;

rho = abs(U(:, :, 1));
t1  = U(:, :, 2);
t2  = U(:, :, 3);
m   = U(:, :, 4);


u   = t1./rho;
v   = t2./rho;
p   = (gamma-1)*(m-0.5.*rho.*u.^2-0.5.*rho.*v.^2);

% rho = 0.5*(abs(rho)+rho);
% p  = 0.5*(abs(p)+p);

W(:, :, 1)   = rho;
W(:, :, 2)   = u;
W(:, :, 3)   = v;
W(:, :, 4)   = abs(p);

end
