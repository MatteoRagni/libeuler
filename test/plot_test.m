clear all
close all
clc

%% Testing Model

sim test_eulerm

%% Plotting Results

euleri = load('../euleri.csv', '-ascii');
eulere = load('../eulere.csv', '-ascii');

figure('color', 'w');
subplot(2,1,1);
  plot(eulerm(:,1), eulerm(:,3), 'g-', 'displayname', 'simulink');
  hold on;
  plot(eulere(:,1), eulere(:,3), 'r-', 'displayname', 'explicit');
  plot(euleri(:,1), euleri(:,3), 'k-', 'displayname', 'implicit');
  ylabel('x_1');
  yyaxis right;
  plot(eulere(:,1), eulere(:,2), '--', 'displayname', 'u');
  ylabel('u');
  ylim([4,11]);
  grid on;
  
  title('$\dot{x}_1(t) = a_1 (a_2 u(t) - a_3 \sqrt{\gamma x_1(t)}),\quad \dot{x}_2(t) = b_1 (b_2 \sqrt{\gamma x_1(t)} - b_3 \sqrt{\gamma x_2(t)}),\quad T_s = 10^{-2}$', 'interpreter', 'latex')
subplot(2,1,2);
  plot(eulere(:,1), eulere(:,4), 'r-', 'displayname', 'explicit');
  hold on;
  plot(euleri(:,1), euleri(:,4), 'k-', 'displayname', 'implicit');
  ylabel('x_2');
  yyaxis right;
  plot(eulere(:,1), eulere(:,2), '--', 'displayname', 'u');
  legend();
  ylabel('u');
  xlabel('t');
  ylim([4,11]);
  grid on;
  
figure('color', 'w');
subplot(2,1,1);
  plot(euleri(:,1), euleri(:,3) - eulere(:,3), 'k-', 'displayname', 'error (e/i)');
  hold on;
  plot(euleri(:,1), euleri(:,3) - eulerm(:,3), 'b-', 'displayname', 'error (s/i)');
  ylabel('Errors x_1');
  legend();
  yyaxis right;
  plot(eulere(:,1), eulere(:,2), '--', 'displayname', 'u');
  ylabel('u');
  ylim([4,11]);
  grid on;
  
  title('$\dot{x}_1(t) = a_1 (a_2 u(t) - a_3 \sqrt{\gamma x_1(t)}),\quad \dot{x}_2(t) = b_1 (b_2 \sqrt{\gamma x_1(t)} - b_3 \sqrt{\gamma x_2(t)}),\quad T_s = 10^{-2}$', 'interpreter', 'latex')
subplot(2,1,2);
  plot(euleri(:,1), euleri(:,4) - eulere(:,4), 'k-', 'displayname', 'error (e/i)');
  hold on;
  plot(euleri(:,1), euleri(:,4) - eulerm(:,4), 'b-', 'displayname', 'error (s/i)');
  ylabel('Errors x_2');
  legend();
  yyaxis right;
  plot(eulere(:,1), eulere(:,2), '--', 'displayname', 'u');
  ylabel('u');
  ylim([4,11]);
  grid on;
  
  