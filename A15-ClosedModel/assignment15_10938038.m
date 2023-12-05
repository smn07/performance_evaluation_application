% Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

k = 6; %number of stations
N = 80;
Z = 40; %think time

% Probabilites
p12 = 1;
p21 = 0.1;
p23 = 0.4;
p24 = 0.5;
p35 = 0.6;
p36 = 0.4;
p52 = 1;
p62 = 1;
p42 = 1;

% Average service time
AST = [40,50/1000,2/1000,80/1000,80/1000,120/1000];

% P matrix
P = [0, p12, 0, 0, 0, 0;
     0, 0, p23, p24, 0, 0;
     0, 0, 0, 0, p35, p36;
     0, p42, 0, 0, 0, 0;
     0, p52, 0, 0, 0, 0;
     0, p62, 0, 0, 0, 0];

% Vector l for calculating vk (by using v = l*(I-P)^-1)
l = [1,0,0,0,0,0];

% Visits to the stations
v = l * inv(eye(6) - P);

% Demand of the stations
D = v .* AST;

% Throughput
Qk = zeros(k);
for n = 1:N
    for j = 1:k
       if j == 1
           R(j) = 0;
       else
           R(j) = D(j) * (1 + Qk(j));
       end
    end
    Thr = n / (Z + sum(R));
    Qk = Thr .* R;
end

% Average response time
ART = sum(R);

% Utilization vector (for each station)
U = Thr .* D;

% PRINT VALUES
disp("Throughput: ");
disp(Thr);

disp("Average response time: ");
disp(ART);

disp("Utilization vector: ");
disp(U);










