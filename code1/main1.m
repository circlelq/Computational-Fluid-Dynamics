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
x          = linspace(0,1,N)'; % grid
dx         = (x(end)-x(1))/(N-1);% grid spacing
x_dis      = 0.3; % discontinuity surface location
n_dis      = x_dis*N;

U_l  = [1, 0, 2.5];
U_r  = [0.125, 0, 0.25];
U    = zeros(N,3);

for i=1:N
    if i<=n_dis
        U(i,:) = U_l;
    else
        U(i,:) = U_r;
    end
end
U2           = U;
U21          = U;
W            = U2W(U); % 原始变量 rho,u,p
F            = zeros(N,3); % flux vector​
Fhat         = W2F(W);
steps        = 0;
flag         = 1;
current_time = 0;
t_max        = 0.2;
maxSteps     = 3e2;


aviobj=VideoWriter('2Mac','MPEG-4');
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
    MacCormach();

    % Roe scheme
    % Roe();

    
    W = U2W(U);
    
    % visualization
    subplot(3,1,1);
    plot(x,W(:,1),'bo-')
    xlabel('$x$','interpreter','latex');ylabel('density $\rho$','interpreter','latex');
    title(sprintf('MacCormach scheme\n$x$ discontinue = %0.2f, time = %.3f, CFL = %.2f, $N$ = %d',x_dis,current_time,CFL,N),'interpreter','latex')
    subplot(3,1,2);
    plot(x,W(:,2),'bo-');
    xlabel('$x$','interpreter','latex');ylabel('velosity $u$','interpreter','latex');
    subplot(3,1,3);
    plot(x,W(:,3),'bo-');
    xlabel('$x$','interpreter','latex');ylabel('pressure $p$','interpreter','latex');

    drawnow
    currFrame = getframe(fig);
    writeVideo(aviobj,currFrame);
end
close(aviobj); %关闭
