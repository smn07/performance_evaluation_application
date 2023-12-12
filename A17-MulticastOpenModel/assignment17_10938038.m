% Simone Di Ienno, matricola:225606, codice ID: 10938038

clear all;
clc;

l1a = 2/60;
l1b = 3/60;
l1c = 2.5/60;
lambda = [l1a,l1b,l1c];
lambdaSum = l1a + l1b + l1c;

% minutes
S1C = [10,4,6];
S2C = [12,3,6];

% Visits (for each class)
v1c = [l1a/lambdaSum,l1b/lambdaSum,l1c/lambdaSum];
v2c = v1c;

% Throughput
X1c = v1c * lambdaSum;
X2c = v2c * lambdaSum;

% Utilization
U1 = sum(X1c .* S1C);
U2 = sum(X2c .* S2C);

% Avoiding to calulate P matrix
D1 = 1 .* S1C; %%v = 1
D2 = 1 .* S2C;

% Average response time per product type
R1c = D1 / (1-U1);
R2c = D2 / (1-U2);
Rc = R1c + R2c;

% Average number of jobs per product type
Nc = lambda .* Rc;

% Class-independent average system response time
R = sum((lambda/lambdaSum) .* Rc);


% PRINT VALUES
disp("Utilizations: ");
disp(U1);
disp(U2);
disp("Average response time per product type: ");
disp(Rc);
disp("Average number of jobs per product type: ");
disp(Nc);
disp("Class-independent average system response time: ");
disp(R);

