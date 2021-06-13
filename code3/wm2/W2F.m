function [ F ] = W2F( W )
% transform W to F

global gamma;

rho  = W(:, :, 1);
u    = W(:, :, 2);
v    = W(:, :, 3);
p    = W(:, :, 4);

F(:, :, 1) = rho.*u;
F(:, :, 2) = rho.*u.^2+p;
F(:, :, 3) = rho.*u.*v;
F(:, :, 4) =  p.*u*gamma/(gamma-1)+0.5*rho.*u.*(u.^2+v.^2);

end
