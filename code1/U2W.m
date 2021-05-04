function [ W ] = U2W( U )
% transform U to W
%   U: [rho rho*u rho*E]
%   W: [rho u p]

global gamma;

W   = U(:,1);
t   = U(:,2);
m   = U(:,3);

rho = abs(W);
u   = t./W;
p   = abs((gamma-1)*(m-0.5.*rho.*u.^2));

% rho = 0.5*(abs(rho)+rho);
% p  = 0.5*(abs(p)+p);

W   = [rho,u,p];

end