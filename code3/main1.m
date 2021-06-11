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
% plot(x/D,U0_axi,'-');
% xlabel('$x/D$','interpreter','latex');
% ylabel('$U_0$','interpreter','latex');

% saveas(gcf,'Re_80_U0_axi','epsc')

global gamma
global dt
global dx dy
global Nx Ny
global U
global F
global W
global Fhat Ghat
global U2
global U21
gamma      = 1.4;
% R        = 287;
Nx         = 300; % x grid number
Ny         = 100; % y grid number
CFL        = 0.8; % CFL number
x          = linspace(0, 3, Nx)'; % x grid
dx         = (x(end)-x(1))/(Nx-1);% x grid spacing
y          = linspace(0, 1, Ny)'; % y grid
dy         = (y(end)-y(1))/(Ny-1);% y grid spacing

U0         = [1.4, 3, 0, 1];
U          = zeros(Nx, Ny, 4);

for i = 1:Nx
    for j = 1:Ny
        U(i,j,:) = U0;
    end
end

U2           = U;
U21          = U;
W            = U2W(U); % 原始变量 rho, u, v, p
F            = zeros(Nx, Ny, 4); % flux vector
Fhat         = W2F(W);
Ghat         = W2G(W);
steps        = 0;
flag         = 1;
current_time = 0;
t_max        = 4;
maxSteps     = 3e2;
% U2=W2U(W);


% aviobj=VideoWriter('Roe','MPEG-4');
% open(aviobj);
fig = figure;
set(gcf,'unit','centimeters','position',[20 20 90 30])

while (flag)
    steps   = steps+1;
    dt = CFL*dx/max(max(abs(W(:, :, 2))+sqrt(gamma*W(:, :, 3)./W(:, :, 1))));
    % timestep
    
    if steps > maxSteps % stoping criteria
        disp('WARNING: maxSteps reached!');
        break;
    end
    
    if current_time+dt >= t_max % stoping criteria
        dt = t_max-current_time;
        flag = 0;
    end
    current_time = current_time + dt;
    F = W2F(W);
    G = W2G(W);

    % Roe scheme
    % Roe();

    
    W = U2W(U);
    
    % visualization
    % subplot(3,1,1);
    contour(W(:, :, 1))
    xlabel('$x$','interpreter','latex');
    ylabel('density $\rho$','interpreter','latex');
    % title(sprintf('Roe scheme\n$x$ discontinue = %0.2f, time = %.3f, CFL = %.2f, $Nx$ = %d',x_dis,current_time,CFL,Nx),'interpreter','latex')
    % subplot(3,1,2);
    % plot(x,W(:, :, 2),'bo-');
    % xlabel('$x$','interpreter','latex');ylabel('velosity $u$','interpreter','latex');
    % subplot(3,1,3);
    % plot(x,W(:, :, 3),'bo-');
    % xlabel('$x$','interpreter','latex');ylabel('pressure $p$','interpreter','latex');
    view(90,-90)
    drawnow
    % currFrame = getframe(fig);
    % writeVideo(aviobj,currFrame);
end
% close(aviobj); %关闭
