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
gamma    = 1.4;
R        = 287;

N          = 300; % grid number
CFL        = 0.1; % CFL number
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
W            = U2W(U);
F            = zeros(N,3); % flux vector​
steps        = 0;
flag         = 1;
current_time = 0;
t_max        = 0.2;
maxSteps     = 1e4;


aviobj=VideoWriter('1MacC.avi');
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
    % for i = 2:N-1

    %     % U2(i,:) = 0.5*(U(i+1,:) + U(i-1,:)) - 0.5*dt/dx * (F(i+1,:) - F(i-1,:));

    % end
    
    % MacCormach scheme
    for i = 2:N-1

        U21(i,:) = U(i,:) - dt/dx * (F(i,:)-F(i-1,:));
    end
    F = W2F(U2W(U21));
	for i = 2:N-1
        U2(i,:) = 0.5 * (U(i,:)+U21(i,:)) - 0.5 * dt/dx * (F(i+1,:)-F(i,:));
    end

    U = U2;
    
    W = U2W(U);
    
    % visualization
    subplot(3,1,1);
    plot(x,W(:,1),'bo-')
    xlabel('$x$','interpreter','latex');ylabel('density $\rho$','interpreter','latex');
    title(sprintf('$x$ discontinue = %0.2f, time = %.3f, CFL = %.2f, $N$ = %d',x_dis,current_time,CFL,N),'interpreter','latex')
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


function [ W ] = U2W( U )
% transform U to W
%   U: [rho rho*u rho*E]
%   W: [rho u p]

global gamma;

W   = U(:,1);
t   = U(:,2);
m   = U(:,3);

rho = W;
u   = t./W;
p   = (gamma-1)*(m-0.5.*rho.*u.^2);

% rho = 0.5*(abs(rho)+rho);
% p  = 0.5*(abs(p)+p);

W   = [rho,u,p];

end



function [ U ] = W2U( W )
% tranWform W to dU
%   W: [rho u p]
%   U: [rho rho*u rho*E]

global gamma;


rho = W(:,1);
u   = W(:,2);
p   = W(:,3);

U  = [rho,rho.*u,p/(gamma-1)+0.5*rho.*u.^2];

end

function [ F ] = W2F( W )
% transform W to F
%   W: [rho u p]
%   F: [rho*u rho*u^2+p rho*u*H]

global gamma;

rho=W(:,1);
u=W(:,2);
p=W(:,3);

F=[rho.*u,rho.*u.^2+p,p.*u*gamma/(gamma-1)+0.5*rho.*u.^3];

end

