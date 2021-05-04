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
global dx
global N
global U
global F
global W
global Fhat
global U2
global U21
gamma      = 1.4;
% R        = 287;
N          = 1000; % grid number
CFL        = 0.8; % CFL number
x          = linspace(-5,5,N)'; % grid
dx         = (x(end)-x(1))/(N-1);% grid spacing
x_dis      = -4; % discontinuity surface location
n_dis      = N/10;

W_l  = [3.857143, 2.629369, 10.33333];
W    = zeros(N,3);

for i=1:N
    if i<=n_dis
        W(i,:) = W_l;
    else
        W(i,1) = 1 + 0.2*sin(5*x(i));
        W(i,2) = 0;
        W(i,3) = 1;
    end
end
U            = W2U(W); % 原始变量 rho,u,p
F            = zeros(N,3); % flux vector​
Fhat         = W2F(W);
steps        = 0;
flag         = 1;
current_time = 0;
t_max        = 1.8;
maxSteps     = 1e4;
U2           = U;
U21          = U;

aviobj=VideoWriter('4Roe','MPEG-4');
open(aviobj);
fig = figure;
set(gcf,'unit','centimeters','position',[20 20 30 30])

while (flag)
    steps   = steps+1;
    dt = CFL*dx/max(abs(W(:,2))+sqrt(gamma*W(:,3)./W(:,1)));
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

    % LF scheme
    % LF();
    
    % MacCormach scheme
    % MacCormach();

    % Roe scheme
    Roe();

    % boundary conditions
    U(2,:) = U(1,:);
    U(end-1,:) = U(end,:);
    
    W = U2W(U);
    
    % visualization
    subplot(3,1,1);
    plot(x,W(:,1),'bo-')
    xlabel('$x$','interpreter','latex');ylabel('density $\rho$','interpreter','latex');
    xticks(-5:1:5)
    yticks(0:2:6)
    ylim([0 6])
    title(sprintf('Roe scheme\n$x$ discontinue = %0.2f, time = %.3f, CFL = %.2f, $N$ = %d',x_dis,current_time,CFL,N),'interpreter','latex')
    subplot(3,1,2);
    plot(x,W(:,2),'bo-');
    xlabel('$x$','interpreter','latex');ylabel('velosity $u$','interpreter','latex');
    xticks(-5:1:5)
    ylim([0 4])
    yticks(0:1:4)
    subplot(3,1,3);
    plot(x,W(:,3),'bo-');
    xlabel('$x$','interpreter','latex');ylabel('pressure $p$','interpreter','latex');
    xticks(-5:1:5)
    yticks(0:3:15)
    ylim([0 15])

    drawnow
    currFrame = getframe(fig);
    writeVideo(aviobj,currFrame);
end
close(aviobj); %关闭
