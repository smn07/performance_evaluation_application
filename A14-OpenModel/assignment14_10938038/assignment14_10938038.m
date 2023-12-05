% Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

% Arrival rates
l2 = 2;
l1 = 3;
lambda = l1 + l2;

% Probabalities considering each station
p12 = 0.8;
p20 = 0.2;
p23 = 0.3;
p24 = 0.5;
p32 = 1;
p42 = 1;

% P matrix
P = [0, p12, 0, 0;
     0, 0, p23, p24;
     0, p32, 0, 0;
     0, p42, 0, 0];

% l vector -> external arrival from one station divided by the total of
% arrivals
l = [3/lambda, 2/lambda, 0, 0];

% Visits to the three stations
v = l * inv(eye(4) - P);

% Average Service Time
AST = [2, 30/1000, 100/1000, 80/1000];

% Demand of the four stations
D = v .* AST;

% Throughput
Thr = lambda;

% Utilizations
Uk = [lambda*D(1), lambda*D(2), lambda*D(3), lambda*D(4)];

% System response time
Rk = [D(1), D(2)/(1-Uk(2)), D(3)/(1-Uk(3)), D(4)/(1-Uk(4))];
R = sum(Rk);

% Average number of jobs
N = lambda * R;

% PRINT VALUES
disp("Throughput: ");
disp(Thr);
disp("System response time: ");
disp(R);
disp("Average number of jobs: ");
disp(N);




