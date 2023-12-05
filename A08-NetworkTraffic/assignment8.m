%Simone Di Ienno, matricola: 225606, codice ID: 10938038
clear all;
clc;


l1 = 0.33; %low->medium
l2 = 0.4; %medium->high
l3 = 0.05; %high->down
l4 = 0.05; %medium->down
l5 = 0.05; %low->down

u1 = 0.6; %medium->low
u2 = 1; %high->medium
u3 = 6; %down->low
u4 = 6; %down->medium
u5 = 6; %down->high

p1 = 0.6; %down->low
p2 = 0.3; %down->medium
p3 = 0.1; %down->high
%infinitesimal generator
Q = [-l1-l5, l1 , 0 , l5;
    u1 ,-u1-l2-l4, l2, l4;
    0 , u2 ,-u2-l3, l3;
    p1*u3 , p2*u4 , p3*u5 ,-p1*u3-p2*u4-p3*u5];

disp(Q);

% Show on a plot, the evolution of the states of the system starting from the MEDIUM traffic state, in time interval T = [0, 8].
initial1 =[0,1,0,0];%starting from medium (2)
[t, Sol1]=ode45(@(t,x) Q'*x, [0 8], initial1);
figure;
plot(t, Sol1, "-");
title('Plot starting from Medium');
legend('low','medium','high','down')
set(gcf, 'Name', 'Starting from MEDIUM');

% Show on a plot, the evolution of the states of the system starting from the DOWN state, in time interval T = [0, 8].
initial2=[0,0,0,1];%starting from down (4)
[t, Sol2]=ode45(@(t,x) Q'*x, [0 8], initial2);
figure;
plot(t, Sol2, "-");
title('Plot starting from Down');
legend('low','medium','high','down')
set(gcf, 'Name', 'starting from DOWN');




