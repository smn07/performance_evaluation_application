% Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

% FIRST CASE
l1 = 10;
mu1 = 50;
mu2 = 5;
p1 = 0.8;

D = p1 / mu1 + (1-p1) / mu2;
rho1 = l1 * D;
m2 = 2 * (p1 / mu1^2 + (1-p1) / mu2^2);

% Utilization of the system
U1 = l1 * D;

% Average response time
ART1 = D + l1*m2 / 2 / (1-rho1);

% Average number of jobs
ANJ1 = rho1 + l1^2*m2 / 2 / (1-rho1);



% SECOND CASE
l2 = 240;
k = 5;
c = 3;

L2 = l2 / k;
rho2 = D * L2 / c;
ca = 1 / sqrt(k);
sigma = 2 * (p1 / mu1^2 + (1-p1) / mu2^2) - (p1 / mu1 + (1-p1) / mu2)^2;
cv = sqrt(sigma) / D;

sum = 0;
for i = 0 : c-1
    sum = sum + ((c * rho2)^i) / factorial(i);
end
expTheta = (D / (c * (1-rho2))) / (1 + (1-rho2) * (factorial(c) / (c * rho2)^c) * sum);

U2 = L2 * D;

% Average utilization of the system
average_U2 = U2 / c;

% Average response time
ART2 = D + ((ca^2 + cv^2)/2) * expTheta;

% Average number of jobs
ANJ2 = L2 * ART2;



% PRINT VALUES
disp("==== SCENARIO 1 ====")
disp("Utilization of the system: ")
disp(U1);
disp("Average response time: ")
disp(ART1);
disp("Average number of jobs: ")
disp(ANJ1);

disp("==== SCENARIO 2 ====")
disp("Average utilization of the system: ")
disp(average_U2);
disp("Average response time: ")
disp(ART2);
disp("Average number of jobs: ")
disp(ANJ2);





