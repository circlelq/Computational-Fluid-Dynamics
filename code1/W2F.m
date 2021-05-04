function [ F ] = W2F( W )
% transform W to F
%   W: [rho u p]
%   F: [rho*u rho*u^2+p rho*u*H]

global gamma;

rho  = (W(:,1));
u    = W(:,2);
p    = (W(:,3));

F=[rho.*u,rho.*u.^2+p,p.*u*gamma/(gamma-1)+0.5*rho.*u.^3];

end