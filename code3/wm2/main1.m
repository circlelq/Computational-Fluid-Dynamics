%
%  main1.m
%  code3
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/6/9.
%

clear;clc;close all;

global gamma
global dx dy
global Nx Ny
global F G
global Fhat Ghat
global CFL
gamma      = 1.4;

Ny         = 250; % y grid number
Nx         = Ny*3; % x grid number
CFL        = 0.8; % CFL number
x          = linspace(0, 3, Nx)'; % x grid
dx         = (x(end)-x(1))/(Nx-1);% x grid spacing
y          = linspace(0, 1, Ny)'; % y grid
dy         = (y(end)-y(1))/(Ny-1);% y grid spacing
load('W/3.7988W.mat')
% W0         = [1.4, 3, 0, 1];
% W          = zeros(Nx, Ny, 4);

% for i = 1:Nx
%     for j = 1:Ny
%         W(i,j,:) = W0;
%     end
% end

dt           = 0;
dT           = 4; % every * time to record
U            = W2U(W);
F            = zeros(Nx, Ny, 4); % x flux vector
G            = zeros(Nx, Ny, 4); % y flux vector
Fhat         = W2F(W);
Ghat         = W2G(W);
steps        = 0;
flag         = 1;
current_time = 3.7988;
t_max        = 8;
h = (Ny*0.2);
d = (Nx*0.2);


while (flag)
    steps   = steps+1;

    % timestep
    disp(current_time);
    W = U2W(U);
    
    % boundary conditon
    W = bc(W);

    U = W2U(W);
    
    if current_time+dt >= t_max % stoping criteria
        dt = t_max-current_time;
        disp('Finished!');
        flag = 0;
    end

    current_time = current_time + dt;

    F = W2F(W);
    G = W2G(W);
    

    % Roe scheme
    [tempflux, dt] = Roe(U);
    tempU1 = U - dt * tempflux;
    [tempflux, ~] = Roe(tempU1);
    tempU2 = 0.75 * U + 0.25 * tempU1 - 0.25 * dt * tempflux;
    [tempflux, ~] = Roe(tempU2);
    U = 1/3 * U + 2/3 * tempU2 - 2/3 * dt * tempflux;
    

    % save data

    if abs(current_time - dT) < 2e-3
        save(strcat('W/',num2str(current_time),'W.mat'),'W')
        dT = dT + 0.2;
    end
end       

save W.mat W
