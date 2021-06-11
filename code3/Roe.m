%
%  Roe.m
%  code3
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/6/9.
%

function [] = Roe()
	global dt
	global dx
    global Nx Ny
    global U
	global F
	global W
	global Fhat Ghat
	global gamma

	A = zeros(4,4);
	B = zeros(4,4);
	Lambda_x = zeros(4,4);
	Lambda_y = zeros(4,4);
	R_x = zeros(4,4);
	R_x(1,1) = 1;
	R_x(1,3) = 1;
	R_x(1,4) = 1;
	R_x(3,2) = 1;

	R_y = zeros(4,4);
	R_y(1,2) = 1;
	R_y(1,3) = 1;
	R_y(1,4) = 1;
	R_y(2,1) = 1;

    for i = 1:Nx-1
    	for j = 1:Ny-1
	    	% rhobar = (0.5*sqrtRho_x)^2;
	    	sqrtRho_x = sqrt(W(i, j, 1))+sqrt(W(i+1, j, 1));
	    	ubar_x = (sqrt(W(i, j, 1))*W(i, j, 2)+sqrt(W(i+1, j, 1))*W(i+1, j, 2))/sqrtRho_x;
	    	vbar_x = (sqrt(W(i, j, 1))*W(i, j, 3)+sqrt(W(i+1, j, 1))*W(i+1, j, 3))/sqrtRho_x;
	    	H_x1 = gamma*W(i, j, 4)/((gamma-1)*W(i, j, 1))+0.5*W(i, j, 2)^2+0.5*W(i, j, 3)^2;
	    	H_x2 = gamma*W(i+1, j, 4)/((gamma-1)*W(i+1, j, 1))+0.5*W(i+1, j, 2)^2+0.5*W(i+1, j, 3)^2;
	    	Hbar_x = (sqrt(W(i, j, 1))*H_x1+sqrt(W(i+1, j, 1))*H_x2)/sqrtRho_x;
	    	abar_x = sqrt((gamma-1)*(Hbar_x-0.5*ubar_x^2-0.5*vbar_x^2));
	    	R_x(2,1) = ubar_x;
	    	R_x(2,3) = ubar_x - abar_x;
	    	R_x(2,4) = ubar_x + abar_x;
	    	R_x(3,3) = vbar_x;
	    	R_x(3,4) = vbar_x;
	    	R_x(4,1) = 0.5*ubar_x^2-0.5*vbar_x^2;
	    	R_x(4,2) = vbar_x;
	    	R_x(4,3) = Hbar_x - ubar_x * abar_x;
	    	R_x(4,4) = Hbar_x + ubar_x * abar_x;
	    	Lambda_x(1,1) = abs(ubar_x);
	    	Lambda_x(2,2) = abs(ubar_x);
	    	Lambda_x(3,3) = abs(ubar_x-abar_x);
	    	Lambda_x(4,4) = abs(ubar_x+abar_x);
	    	A = R_x*Lambda_x*inv(R_x);


	    	sqrtRho_y = sqrt(W(i, j, 1))+sqrt(W(i, j+1, 1));
	    	ubar_y = (sqrt(W(i, j, 1))*W(i, j, 2)+sqrt(W(i, j+1, 1))*W(i, j+1, 2))/sqrtRho_y;
	    	vbar_y = (sqrt(W(i, j, 1))*W(i, j, 3)+sqrt(W(i, j+1, 1))*W(i, j+1, 3))/sqrtRho_y;
	    	H_y1 = gamma*W(i, j, 4)/((gamma-1)*W(i, j, 1))+0.5*W(i, j, 2)^2+0.5*W(i, j, 3)^2;
	    	H_y2 = gamma*W(i, j+1, 4)/((gamma-1)*W(i, j+1, 1))+0.5*W(i, j+1, 2)^2+0.5*W(i, j+1, 3)^2;
	    	Hbar_x = (sqrt(W(i, j, 1))*H_y1+sqrt(W(i, j+1, 1))*H_y2)/sqrtRho_y;
	    	abar_x = sqrt((gamma-1)*(Hbar_x-0.5*ubar_y^2-0.5*vbar_y^2));
	    	R_y(2,3) = ubar_y;
	    	R_y(2,4) = ubar_y;
	    	R_y(3,2) = vbar_y;
	    	R_y(3,3) = vbar_y - abar_y;
	    	R_y(3,4) = vbar_y + abar_y;
	    	R_y(4,1) = ubar_y;
	    	R_y(4,2) = 0.5*vbar_y^2-0.5*ubar_y^2;
	    	R_y(4,3) = Hbar_y - vbar_y * abar_y;
	    	R_y(4,4) = Hbar_y + vbar_y * abar_y;
	    	Lambda_y(1,1) = abs(vbar_y);
	    	Lambda_y(2,2) = abs(vbar_y);
	    	Lambda_y(3,3) = abs(vbar_y-abar_y);
	    	Lambda_y(4,4) = abs(vbar_y+abar_y);
	    	B = R_y*Lambda_y*inv(R_y);

	        Fhat(i, j, :) = (F(i, j, :)+F(i+1, j, :))/2 - 0.5 * (A * (U(i+1, j, :)-U(i, j, :))')';
	        Ghat(i, j, :) = (G(i, j, :)+G(i, j+1, :))/2 - 0.5 * (B * (U(i, j+1, :)-U(i, j, :))')';
	    end
    end

	for i = 2:N-1
        U(i, j, :) = U(i, j, :) - dt/dx * (Fhat(i, j, :)-Fhat(i-1,:));
    end

end

