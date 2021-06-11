function [ W ] = U2W( U )
% transform U to W
%   U: [rho rho*u rho*E]
%   W: [rho u v p]

global gamma;

W   = U(:, :, 1);
t1  = U(:, :, 2);
t2  = U(:, :, 3);
m   = U(:, :, 4);


rho = W;
u   = t1./W;
v   = t2./W;
p   = (gamma-1)*(m-0.5.*rho.*u.^2-0.5.*rho.*v.^2);

% rho = 0.5*(abs(rho)+rho);
% p  = 0.5*(abs(p)+p);

W(:, :, 1)   = rho;
W(:, :, 2)   = u;
W(:, :, 3)   = v;
W(:, :, 4)   = p;

end
