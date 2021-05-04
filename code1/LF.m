function [] = LF()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
	global dt
	global dx
    global N
    global U
	global F
	global U2
	global U21
    for i = 2:N-1

        U2(i,:) = 0.5*(U(i+1,:) + U(i-1,:)) - 0.5*dt/dx * (F(i+1,:) - F(i-1,:));

    end
    U = U2;
end

