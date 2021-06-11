function [ U ] = W2U( W )
% tranWform W to U

global gamma;

rho  = W(:, :, 1);
u    = W(:, :, 2);
v    = W(:, :, 3);
p    = W(:, :, 4);

U(:, :, 1) = rho;
U(:, :, 2) = rho.*u;
U(:, :, 3) = rho.*v;
U(:, :, 4) = p/(gamma-1)+0.5*rho.*(u.^2+v.^2);

end
