function [ G ] = W2F( W )
% transform W to G

global gamma;

rho  = W(:, :, 1);
u    = W(:, :, 2);
v    = W(:, :, 3);
p    = W(:, :, 4);

G(:, :, 1) = rho.*v;
G(:, :, 2) = rho.*u.*v;
G(:, :, 3) = rho.*v.^2+p;
G(:, :, 4) = p.*v*gamma/(gamma-1)+0.5*rho.*v.*(u.^2+v.^2);

end
