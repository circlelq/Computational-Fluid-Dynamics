function [ U ] = W2U( W )
% tranWform W to dU
%   W: [rho u p]
%   U: [rho rho*u rho*E]

global gamma;


rho  = (W(:,1));
u    = W(:,2);
p    = (W(:,3));

U  = [rho,rho.*u,p/(gamma-1)+0.5*rho.*u.^2];

end