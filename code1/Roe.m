function [] = Roe()
	global dt
	global dx
    global N
    global U
	global F
	global W
	global Fhat
	global gamma
	A = zeros(3,3);
	Lambda = zeros(3,3);
	R = zeros(3,3);
	R(1,1) = 1;
	R(1,2) = 1;
	R(1,3) = 1;
    for i = 1:N-1
    	% rhobar = (0.5*(sqrt(W(i,1))+sqrt(W(i+1,1))))^2;
    	ubar = (sqrt(W(i,1))*W(i,2)+sqrt(W(i+1,1))*W(i+1,2))/(sqrt(W(i,1))+sqrt(W(i+1,1)));
    	H1 = gamma*W(i,3)/((gamma-1)*W(i,1))+0.5*W(i,2)^2;
    	H2 = gamma*W(i+1,3)/((gamma-1)*W(i+1,1))+0.5*W(i+1,2)^2;
    	Hbar = (sqrt(W(i,1))*H1+sqrt(W(i+1,1))*H2)/(sqrt(W(i,1))+sqrt(W(i+1,1)));
    	abar = sqrt((gamma-1)*(Hbar-0.5*ubar^2));
    	R(2,1) = ubar - abar;
    	R(2,2) = ubar;
    	R(2,3) = ubar + abar;
    	R(3,1) = Hbar - ubar * abar;
    	R(3,2) = 0.5*ubar^2;
    	R(3,3) = Hbar + ubar * abar;
    	Lambda(1,1) = abs(ubar-abar);
    	Lambda(2,2) = abs(ubar);
    	Lambda(3,3) = abs(ubar+abar);
    	A = R*Lambda*inv(R);

        Fhat(i,:) = (F(i,:)+F(i+1,:))/2 - 0.5 * (A * (U(i+1,:)-U(i,:))')';
    end
	for i = 2:N-1
        U(i,:) = U(i,:) - dt/dx * (Fhat(i,:)-Fhat(i-1,:));
    end
end

