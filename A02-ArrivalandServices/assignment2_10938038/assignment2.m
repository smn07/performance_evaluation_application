% SIMONE DI IENNO, MATRICOLA: 225606, CODICE ID: 10938038, ASSIGNMENT2

clear all;

%TRACE 1
trace1 = csvread('Trace1.csv');

% #arrivals
A = size(trace1,1);
% #completions
C = size(trace1,1);

arrivalstime = cumsum(trace1(:, 1)); %array of arrivals time

completions = sum(trace1(1, 1:2));
for i = 2:A
    completions=[completions;max(arrivalstime(i),completions(i-1))+trace1(i,2)];  % i'm building the completions array
end

% i-th response time
ri = completions - arrivalstime;

% Response time
R = sum(ri)/C;

INv = [arrivalstime, ones(length(arrivalstime),1)];
OUTv = [completions, -ones(length(completions),1)];

matrix = [INv; OUTv]; % matrix with ones and -ones

matrix = sortrows(matrix,1); % sort with respect to the first column
matrix(:,3) = cumsum(matrix(:,2)); % cumulative sum
matrix(1:end-1, 4) = matrix(2:end,1) - matrix(1:end-1,1);

matrix(:,5) = (matrix(:,3) > 0) .* matrix(:,4);

% Busy time
B = sum(matrix(:,5));

% Time
T = completions(end) - arrivalstime(1);

% Utilization time
U = B/T;

% i'm considering only when there are no jobs
matrix(:,5) = (matrix(:,3) == 0);

% frequency idle
frequency = sum(matrix(:,5))/T;

% average idle time
AIT = (T-B)/sum(matrix(:,5));




%TRACE 2
trace2 = csvread('Trace2.csv');

% #arrivals
A2 = size(trace2,1);
% #completions
C2 = size(trace2,1);

arrivalstime2 = cumsum(trace2(:, 1)); %array of arrivals time

completions2 = sum(trace2(1, 1:2));
for i = 2:A
    completions2=[completions2;max(arrivalstime2(i),completions2(i-1))+trace2(i,2)];  % i'm building the completions array
end

% i-th response time
ri2 = completions2 - arrivalstime2;

% Response time
R2 = sum(ri2)/C2;

INv2 = [arrivalstime2, ones(length(arrivalstime2),1)];
OUTv2 = [completions2, -ones(length(completions2),1)];

matrix2 = [INv2; OUTv2]; % matrix with ones and -ones

matrix2 = sortrows(matrix2,1); % sort with respect to the first column
matrix2(:,3) = cumsum(matrix2(:,2)); % cumulative sum
matrix2(1:end-1, 4) = matrix2(2:end,1) - matrix2(1:end-1,1);

matrix2(:,5) = (matrix2(:,3) > 0) .* matrix2(:,4);

% Busy time
B2 = sum(matrix2(:,5));

% Time
T2 = completions2(end) - arrivalstime2(1);

% Utilization time
U2 = B2/T2;

% i'm considering only when there are no jobs
matrix2(:,5) = (matrix2(:,3) == 0);

% frequency idle
frequency2 = sum(matrix2(:,5))/T2;

% average idle time
AIT2 = (T2-B2)/sum(matrix2(:,5));



%TRACE 3
trace3 = csvread('Trace3.csv');

% #arrivals
A3 = size(trace3,1);
% #completions
C3 = size(trace3,1);

arrivalstime3 = cumsum(trace3(:, 1)); %array of arrivals time

completions3 = sum(trace3(1, 1:2));
for i = 2:A
    completions3=[completions3;max(arrivalstime3(i),completions3(i-1))+trace3(i,2)];  % i'm building the completions array
end

% i-th response time
ri3 = completions3 - arrivalstime3;

% Response time
R3 = sum(ri3)/C3;

INv3 = [arrivalstime3, ones(length(arrivalstime3),1)];
OUTv3 = [completions3, -ones(length(completions3),1)];

matrix3 = [INv3; OUTv3]; % matrix with ones and -ones

matrix3 = sortrows(matrix3,1); % sort with respect to the first column
matrix3(:,3) = cumsum(matrix3(:,2)); % cumulative sum
matrix3(1:end-1, 4) = matrix3(2:end,1) - matrix3(1:end-1,1);

matrix3(:,5) = (matrix3(:,3) > 0) .* matrix3(:,4);

% Busy time
B3 = sum(matrix3(:,5));

% Time
T3 = completions3(end) - arrivalstime3(1);

% Utilization time
U3 = B3/T3;

% i'm considering only when there are no jobs
matrix3(:,5) = (matrix3(:,3) == 0);

% frequency idle
frequency3 = sum(matrix3(:,5))/T3;

% average idle time
AIT3 = (T3-B3)/sum(matrix3(:,5));



%PRINT VALUES FOR THE ASSIGNMENT 2
disp('TRACE1.CSV VALUES:');
disp('Average Response Time = ');
disp(R);
disp('Utilization = ');
disp(U);
disp('Frequency the system returns to idle = ');
disp(frequency);
disp('Average Idle Time = ');
disp(AIT);

disp('TRACE2.CSV VALUES:');
disp('Average Response Time = ');
disp(R2);
disp('Utilization = ');
disp(U2);
disp('Frequency the system returns to idle = ');
disp(frequency2);
disp('Average Idle Time = ');
disp(AIT2);

disp('TRACE3.CSV VALUES:');
disp('Average Response Time = ');
disp(R3);
disp('Utilization = ');
disp(U3);
disp('Frequency the system returns to idle = ');
disp(frequency3);
disp('Average Idle Time = ');
disp(AIT3);
