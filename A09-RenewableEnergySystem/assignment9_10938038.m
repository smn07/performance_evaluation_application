%Simone Di Ienno; matricola: 225606; codice ID: 10938038
clear all;
clc;


l2 = 1/18; %from state "night" to "scanning"
l3 = 1/3; %from state "sunny day" to "scanning"
l4 = 1/8; %from state "cloudy day" to "scanning"
m = 1/2; % from "scanning" to each state

p1 = 0.5; % from scanning to night
p2 = 0.33; % from scanning to sunny day
p3 = 0.17; % from scanning to cloudy day

Q = [-p1*m-p2*m-p3*m, p1*m, p2*m, p3*m;
        l2, -l2,   0  , 0;
        l3, 0, -l3, 0;
        l4, 0, 0, -l4];

initial_state = [0, p1, p2, p3];%initial state depending on probabilities

rew_power_consump = [12, 0.1, 0.1, 0.1]; % rewards for power consumption
rew_utiliz = [1, 0, 0, 0];% rewards for utilization (i consider only scanning state)
rew_thr=[0,0,0,0; % rewards for throughput -> this is a matrix
         1,0,0,0;
         1,0,0,0;
         1,0,0,0];

[t, Sol1] = ode45(@(t,x) Q'*x, [0 1440], initial_state'); % for 24 hours
avg_power_consump = Sol1(end,:) * rew_power_consump';
utilization = Sol1(end,:) * rew_utiliz';
throughput= Sol1(end,:) * (sum((rew_thr .* Q)'))';

% PRINT VALUES
disp('Q matrix:');
disp(Q);
disp('Reward vector for average power consumption:');
disp(rew_power_consump);
disp('Reward vector for utilization:');
disp(rew_utiliz);
disp('Reward matrix for throughput:');
disp(rew_thr);
disp('Average power consumption:');
disp(avg_power_consump);
disp('Average utilization:');
disp(utilization);
disp('Throughput:');
disp(throughput*1440);
