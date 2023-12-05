% Simone Di Ienno, matricola:225606, codice ID: 10938038

clear all;
clc;

%lambda
l = 150/60; % req/seconds
D = 350/10^3; % seconds

k = 32;
rho = l * D;

% Utilization
U = (rho - rho^(k+1)) / (1 - rho^(k+1));

% Loss probability
pl = (rho^k - rho^(k+1)) / (1-rho^(k+1));

% Average number of jobs in the system
ANJ = (rho/(1-rho)) - (((k+1)*rho^(k+1)) / (1-rho^(k+1)));

% Drop rate
dr = l * pl;  %(because pl = pk)

% Average Response Time
ART = ANJ / (l * (1-pl)); %(because pl = pk)

% Average time spent in the queue (waiting for service)
ATSQ = ART - D;

% PRINTS MM1K
disp("========MM1K========");
disp("Utilization:");
disp(U);
disp("Loss probability:");
disp(pl);
disp("Average number of jobs in the system:");
disp(ANJ);
disp("Drop rate");
disp(dr);
disp("Average Response Time");
disp(ART);
disp("Average time spent in the queue (waiting for service)");
disp(ATSQ);

% MM2K

%lambda
l2 = 250/60; % req/seconds
c = 2;
rho2 = (l2*D) / c;

% Total utilization of the system
tmp = 0;

for i = 0:(c-1)
	tmp = tmp + (c*rho2)^i / factorial(i);
end

p0 = ((c*rho2)^c*(1-rho2^(k-c+1)) / factorial(c) / (1-rho2) + tmp)^-1;

totalU = 0;

for i=1:c
	totalU = totalU + i * p(i,c,p0,rho2);
end

for i = c+1:k
	totalU = totalU + c * p(i,c,p0,rho2);
end

% Average utilization of the system
averageU = totalU / 2;

% Loss probability
pl2 = p(k,c,p0,rho2);

% Average number of jobs in the system
ANJ2 = 0;
for i=1:k
	ANJ2 = ANJ2 + i * p(i,c,p0,rho2);
end

% Drop rate
dr2 = l2 * pl2;

% Average response time
ART2 = ANJ2 / (l2 * (1 - pl2));

% Average time spent in the queue (waiting for service)
ATSQ2 = ART2 - D;

% PRINTS MM2K
disp("========MM2K========");
disp("Utilization:");
disp(U);
disp("Average utilization:");
disp(averageU);
disp("Loss probability:");
disp(pl2);
disp("Average number of jobs in the system:");
disp(ANJ2);
disp("Drop rate");
disp(dr2);
disp("Average Response Time");
disp(ART2);
disp("Average time spent in the queue (waiting for service)");
disp(ATSQ2);

function F = p(i,c,p0,rho_f)
	if i < c
		F = p0 / factorial(i) * (c*rho_f)^i;
	else
		F = p0*c^c * rho_f^i / factorial(c);
	end
end













