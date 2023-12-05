% Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

% FIRST CASE

% Probabalities considering each station
% 1 = Terminals
% 2 = CPU
% 3 = Disk
% 4 = RAM
p12 = 1;
p21 = 0.1;
p23 = 0.3;
p24 = 0.6;
p32 = 0.85;
p34 = 0.15;
p42 = 0.75;
p43 = 0.25;

% Average service time
AST1 = [10,20/1000,10/1000,3/1000];

% Vector l for calculating vk (by using v = l*(I-P)^-1)
l = [1,0,0,0];

% Matrix P
P = [0, p12, 0, 0;
	 0, 0, p23, p24;
	 0, p32, 0, p34;
	 0, p42, p43, 0];

% Visits to the four stations
v1 = l * inv(eye(4) - P);
% Demand of the four stations
D1 = v1 .* AST1;


% SECOND CASE

lambda = 0.3; %(enter the system)

% 1 = CPU
% 2 = Disk
% 3 = RAM

p12 = 0.3;
p13 = 0.6;
p21 = 0.8;
p23 = 0.15;
p31 = 0.75;
p32 = 0.25;

% Average service time
AST2 = [20/1000,10/1000,3/1000];

P2 = [0, p12, p13;
	 p21, 0, p23;
	 p31, p32, 0];
% Vector l for calculating vk (by using v = l*(I-P)^-1)
l2 = [1,0,0];

% Visits to the three stations
v2 = l2 * inv(eye(3) - P2);

% Demand of hte three stations
D2 = v2 .* AST2;

% Throughput of the three stations
Thr = v2 * lambda;



% PRINT VALUES
disp("==== Scenario 1 ====");
disp("The visits to the four stations: ");
disp(v1);
disp("The demand of the four stations");
disp(D1);

disp("==== Scenario 2 ====");
disp("The visits to the three stations: ");
disp(v2);
disp("The demand of the three stations");
disp(D2);
disp("Throughput of the three stations");
disp(Thr);



