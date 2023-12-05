% Simone Di Ienno, matricola: 225606, codice ID: 10938038 

clear all;
clc;

% MM1
l = 40;
D = 16/10^3; %seconds

% (traffic intensity)
rho = l * D;

% Utilization
U = rho;

% Probability of having exactly one job in the system
p1job = (1 - rho) * rho;

% Probability of having less than 10 jobs in the system
pLess10jobs = 1 - rho^10;

% Average queue length (jobs not in service)
AQL = rho^2 / (1 - rho);

% Average response time
ART = D / (1 - rho);

% Probability that the response time is greater than 0.5 s
pRTgreater05 = exp(-0.5 / ART);

% 90 percentile of the response time distribution
perc90 = -log(1 - 90 / 100) * ART;

% PRINTS FOR MM1
disp("========MM1========");
disp("Utilization:");
disp(U);
disp("Probability of having exactly one job in the system:");
disp(p1job);
disp("Probability of having less than 10 jobs in the system:");
disp(pLess10jobs);
disp("Average queue length (jobs not in service):");
disp(AQL);
disp("Average response time:");
disp(ART);
disp("Probability that the response time is greater than 0.5 s:");
disp(pRTgreater05);
disp("90 percentile of the response time distribution");
disp(perc90);


% MM2
l2 = 90;
rho2 = (l2 * D) / 2;

% Total utilization of the system
totalU = l2 * D;

% Average utilization of the system
averageU = totalU / 2;

% Probability of having exactly one job in the system
p1Job2 = 2 * ((1 - rho2) / (1 + rho2)) * rho2;

% Probability of having less than 10 jobs in the system
pLess10jobs2 = (1 - rho2) / (1 + rho2);
for i=1:9
    pLess10jobs2 = pLess10jobs2 + 2 * (1 - rho2) / (1 + rho2) * rho2^i;
end

%. Average queue length (jobs not in service)
n = 2 * rho2 / (1 - rho2^2);
AQL2 = n - totalU;

% Average response time
ART2 = D / (1 - rho2^2);

% PRINTS FOR MM2
disp("========MM2========");
disp("Utilization:");
disp(totalU);
disp("Average Utilization:");
disp(averageU);
disp("Probability of having exactly one job in the system:");
disp(p1Job2);
disp("Probability of having less than 10 jobs in the system:");
disp(pLess10jobs2);
disp("Average queue length (jobs not in service):");
disp(AQL2);
disp("Average response time:");
disp(ART2);







