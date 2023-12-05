% SIMONE DI IENNO, MATRICOLA: 225606, CODICE PERSONA: 10938038, ASSIGNMENT 3

clear all;

t1 = csvread('Trace1.csv');
t2 = csvread('Trace2.csv');
t3 = csvread('Trace3.csv');

matrix = [t1 t2 t3];
N = length(t1);

% Mean Moment
m1 = sum(matrix) / N;

% Second Moment
m2 = sum(matrix .^ 2) / N;

% Third Moment
m3 = sum(matrix .^ 3) / N;

% Fourth Moment
m4 = sum(matrix .^ 4) / N;

% Variance
c2 = sum((matrix - m1) .^ 2) / N; % var(matrix);

% Third certered moment
c3 = sum((matrix - m1) .^ 3) / N;

% Fourth certered moment
c4 = sum((matrix - m1) .^ 4) / N;

% Sigma (Standard Deviation)
sigma = sqrt(c2);  % std(matrix)

% Skewness
s3 = sum(((matrix - m1) ./ sigma) .^ 3) / N;  % skewness(matrix);

% Fourth stardardized moment
s4 = sum(((matrix - m1) ./ sigma) .^ 4) / N;  % Kurtosis(matrix);

% Excess Kurtosis
exK = s4 - 3;

% Coefficient of variation
cv = sigma ./ m1;

% Median
median = median(matrix);

sortedMatrix = sort(matrix);

h25 = (N-1)*0.25+1;
%integer version of h25:
ih25 = floor(h25);
d25 = h25 - ih25;

h75 = (N-1)*0.75+1;
ih75 = floor(h75);
d75 = h75 - ih75;

% First and third quartile
q1 = sortedMatrix(ih25, :) + (sortedMatrix(ih25+1, :) - sortedMatrix(ih25, :))*d25;
q3 = sortedMatrix(ih75, :) + (sortedMatrix(ih75+1, :) - sortedMatrix(ih75, :))*d75;

% Figure with the Pearson Correlation Coefficient for lags m=1 to m=100
Si = zeros(100, 3);
for i = 1:100
        Si(i, :) = sum((matrix(1:end-i,:)-m1) .* (matrix((i+1):end,:)-m1)) / (N-i);
end

PearsonC = Si ./ c2;

figure;
plot([1:100],PearsonC(:,1),"-")
title('Figure with the Pearson Correlation Coefficient TRACE 1');
grid on;
xlabel('Lag(m)');
ylabel('Pearson Correlation Coefficient');
l = legend;

figure;
plot([1:100],PearsonC(:,2),"-")
title('Figure with the Pearson Correlation Coefficient TRACE 2');
grid on;
xlabel('Lag(m)');
ylabel('Pearson Correlation Coefficient');
l = legend;

figure;
plot([1:100],PearsonC(:,3),"-")
title('Figure with the Pearson Correlation Coefficient TRACE 3');
grid on;
xlabel('Lag(m)');
ylabel('Pearson Correlation Coefficient');
l = legend;

% Approximated CDF of the corresponding distribution
figure;
plot(sortedMatrix(:,1), [1:N]/N, ".")
title('Approximated CDF of the corresponding distribution TRACE 1');
grid on;
xlabel('Lag(m)');
ylabel('Cumulative Probability');
l = legend;

figure;
plot(sortedMatrix(:,2), [1:N]/N, ".")
title('Approximated CDF of the corresponding distribution TRACE 2');
grid on;
xlabel('Lag(m)');
ylabel('Cumulative Probability');
l = legend;

figure;
plot(sortedMatrix(:,3), [1:N]/N, ".")
title('Approximated CDF of the corresponding distribution TRACE 3');
grid on;
xlabel('Lag(m)');
ylabel('Cumulative Probability');
l = legend;



% PRINT VALUES
disp('Mean moment = ');
disp(m1);
disp('Second moment = ');
disp(m2);
disp('Third moment = ');
disp(m3);
disp('Fourth moment = ');
disp(m4);
disp('Variance = ');
disp(c2);
disp('Third centered moment = ');
disp(c3);
disp('Fourth centered moment = ');
disp(c4);
disp('Standard deviation =');
disp(sigma);
disp('Skewness = ');
disp(s3);
disp('Fourth standardized moment = ');
disp(s4);
disp('Excess Kurtosis = ');
disp(exK);
disp('Coefficient of variation = ');
disp(cv);
disp('Median = ');
disp(median);
disp('First quartile = ');
disp(q1);
disp('Third quartile = ');
disp(q3);



