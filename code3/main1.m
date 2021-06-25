%
%  main1.m
%  code3
%
%  Created by Yuan Leiqi (袁磊祺) on 2021/6/9.
%

clear;clc;close all;
set(0,'defaultlinelinewidth',3)
set(0,'defaultaxeslinewidth',2);
set(0,'defaultaxesfontsize',28);
set(0,'defaulttextfontsize',28);
set(0,'DefaultLineMarkerSize',2);
set(0,'Defaultaxesfontname','Times New Roman');
% addpath('./circleData');
% set(gcf,'unit','centimeters','position',[20 20 20 20])
% figure(); box on;
% plot(x/D,W0_axi,'-');
% xlabel('$x/D$','interpreter','latex');
% ylabel('$U_0$','interpreter','latex');

% saveas(gcf,'Re_80_U0_axi','epsc')

global gamma
global dx dy
global Nx Ny
global F G
global Fhat Ghat
global CFL
gamma      = 1.4;

Ny         = 50; % y grid number
Nx         = Ny*3; % x grid number
CFL        = 0.5; % CFL number
x          = linspace(0, 3, Nx)'; % x grid
dx         = (x(end)-x(1))/(Nx-1);% x grid spacing
y          = linspace(0, 1, Ny)'; % y grid
dy         = (y(end)-y(1))/(Ny-1);% y grid spacing

W0         = [1.4, 3, 0, 1];
W          = zeros(Nx, Ny, 4);

for i = 1:Nx
    for j = 1:Ny
        W(i,j,:) = W0;
    end
end

dt           = 0;
dT           = 0; % every 0.1 time to record
U            = W2U(W);
F            = zeros(Nx, Ny, 4); % x flux vector
G            = zeros(Nx, Ny, 4); % y flux vector
Fhat         = W2F(W);
Ghat         = W2G(W);
steps        = 0;
flag         = 1;
current_time = 0;
t_max        = 4;
% maxSteps     = 3e4;
h = (Ny*0.2);
d = (Nx*0.2);

aviobj=VideoWriter('Roe','MPEG-4');
open(aviobj);
fig = figure;
set(gcf,'unit','centimeters','position',[20 20 90 30])

while (flag)
    steps   = steps+1;

    % timestep
    % disp(current_time);
    W = U2W(U);
    
    % boundary conditon
    W = bc(W);

    U = W2U(W);

    % if steps > maxSteps % stoping criteria
    %     disp('WARNING: maxSteps reached!');
    %     break;
    % end
    
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
    

    % visualization

    for i = 1:d
        for j = 2:Ny
            W1(i,j,:) = W(i,j,:);
        end
    end

    for i = d+1:Nx
        for j = h+1:Ny
            W1(i,j,:) = W(i,j,:);
        end
    end
    contourf(W1(:, :, 1))
    xlabel('$y$','interpreter','latex');
    ylabel('$x$','interpreter','latex');
    title(sprintf('Roe scheme, time = %.3f',current_time),'interpreter','latex')

    view(90,-90)
    drawnow
    currFrame = getframe(fig);
    % if abs(current_time - dT) < 1e-2
        writeVideo(aviobj,currFrame);
        % save(strcat(num2str(current_time),'W.mat'),'W')
        % dT = dT + 0.1;
    % end
end
close(aviobj); %关闭
       

% save W.mat W
