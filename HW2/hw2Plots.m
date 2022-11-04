clc; clear all; close all;

%% n=51, dx=0.04 (dt varyinig, dx constant)
dxdt = [0.05 0.1 0.25 0.5 1 1.5 5 10 20];
dt = [0.8 0.4 0.16 0.08 0.04 0.0266667 0.008 0.004 0.002];
err = [0.819062 0.133695 0.0662523 0.0821121 0.0926881 0.0953226 0.100007 0.101538 0.11027];

figure()
subplot(1,2,1)
plot(dxdt, err, '*-','LineWidth',1.5)
hold on 
xline(1, 'r','LineWidth',1.5)
legend('Error vs \Delta t', '\Delta x = \Delta t','FontSize',18)
ylabel('Error','FontSize',18)
xlabel('\Delta x/\Delta t','FontSize',18)
subplot(1,2,2)
plot(dt, err, '*-','LineWidth',1.5)
hold on 
xline(0.04, 'r','LineWidth',1.5)
legend('Error vs \Delta t', '\Delta x = \Delta t','FontSize',18)
ylabel('Error','FontSize',18)
xlabel('\Delta t','FontSize',18)
sgtitle('Error with n=51, \Delta x=0.04','FontSize',24)

figure()
loglog(dt, err, '*-', 'LineWidth',1.5)
hold on
loglog(dt,arrayfun(@(dx) dx^0,dt),'g:','LineWidth',1.5)
loglog(dt,arrayfun(@(dx) dx^1,dt),'r:','LineWidth',1.5)
loglog(dt,arrayfun(@(dx) dx^2,dt),'m:','LineWidth',1.5)
xlabel('dt','FontSize',18)
ylabel('error','FontSize',18)
legend('SL', 'slope 0', 'slope -1', 'slope -2', 'FontSize',24,'Location','southeast')
title('SL Error and Convergence, n=51, \Delta x=0.04','FontSize',24)


%% n=21, dx=0.1 (dt varyinig, dx constant)
% Decided not to use these results as redundant, commenting out
% dxdt = [0.1 0.25 0.5 1 1.5 5 10 20];
% dt=[1 0.4 0.2 0.1 0.0666667 0.02 0.01 0.005];
% err=[2.00317 0.132414 0.126875 0.156255 0.168773 0.182862 0.185717 0.187163];
% 
% figure()
% subplot(1,2,1)
% plot(dxdt, err, '*-','LineWidth',1.5)
% hold on 
% xline(1, 'r','LineWidth',1.5)
% legend('Error vs \Delta t', '\Delta x = \Delta t','FontSize',18)
% ylabel('Error','FontSize',18)
% xlabel('\Delta x/\Delta t','FontSize',18)
% subplot(1,2,2)
% plot(dt, err, '*-','LineWidth',1.5)
% hold on 
% xline(0.1, 'r','LineWidth',1.5)
% legend('Error vs \Delta t', '\Delta x = \Delta t','FontSize',18)
% ylabel('Error','FontSize',18)
% xlabel('\Delta t','FontSize',18)
% sgtitle('Error with n=21, \Delta x=0.1','FontSize',24)
% 
% figure()
% loglog(dt, err, '*-', 'LineWidth',1.5)
% hold on
% loglog(dt,arrayfun(@(dx) dx^0,dt),'g:','LineWidth',1.5)
% loglog(dt,arrayfun(@(dx) dx^1,dt),'r:','LineWidth',1.5)
% loglog(dt,arrayfun(@(dx) dx^2,dt),'m:','LineWidth',1.5)
% xlabel('dt','FontSize',18)
% ylabel('error','FontSize',18)
% legend('SL', 'slope 0', 'slope -1', 'slope -2', Location='northwest')
% title('SL Error and Convergence, n=21, \Delta x=0.1','FontSize',24)

%% dt =0.1 (dt constant, dx varying)
dx=[0.2 0.1 0.0666667 0.05 0.04 0.02 0.01];
dxdt=[2 1 0.666667 0.5 0.4 0.2 0.1];
err=[0.352739 0.156255 0.122079 0.105523 0.0756888 0.0317766 0.0215366];

figure()
subplot(1,2,1)
plot(dxdt, err, '*-','LineWidth',1.5)
hold on 
xline(1, 'r','LineWidth',1.5)
legend('Error vs \Delta x', '\Delta x = \Delta t','FontSize',18)
ylabel('Error','FontSize',18)
xlabel('\Delta x/\Delta t','FontSize',18)
subplot(1,2,2)
plot(dx, err, '*-','LineWidth',1.5)
hold on 
xline(0.1, 'r','LineWidth',1.5)
legend('Error vs \Delta t', '\Delta x = \Delta t','FontSize',18)
ylabel('Error','FontSize',18)
xlabel('\Delta x','FontSize',18)
sgtitle('Error with \Delta t=0.1','FontSize',24)

figure()
loglog(dx, err, '*-', 'LineWidth',1.5)
hold on
loglog(dx,arrayfun(@(dxs) dxs^0,dx),'g:','LineWidth',1.5)
loglog(dx,arrayfun(@(dxs) dxs^1,dx),'r:','LineWidth',1.5)
loglog(dx,arrayfun(@(dxs) dxs^2,dx),'m:','LineWidth',1.5)
xlabel('dx','FontSize',18)
ylabel('error','FontSize',18)
legend('SL', 'slope 0', 'slope -1', 'slope -2', 'FontSize',24,'Location','southeast')
title('SL Error and Convergence \Delta t=0.1','FontSize',24)