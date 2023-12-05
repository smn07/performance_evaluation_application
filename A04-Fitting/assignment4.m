% Simone Di Ienno, matricola: 225606, codice ID: 10938038
clear all;
clc;
t1 = csvread('Trace1.csv');
t2 = csvread('Trace2.csv');

matrix = [t1 t2];
N = length(t1);

sortedMatrix = sort(matrix);

Mu = mean(sortedMatrix);
sigma = std(sortedMatrix);
cv = sigma ./ Mu;

%global variables for weibull and pareto equation functions
global m1_perweibull;
global m2_perweibull;
global m1_perpareto;
global m2_perpareto;

% Mean Moment
m1 = sum(sortedMatrix) / N

% Second Moment
m2 = sum(sortedMatrix .^ 2) / N

% Third Moment
m3 = sum(sortedMatrix .^ 3) / N

% Uniform
disp("Uniform parameters:");
a = m1 - sqrt(12* (m2 - m1.^2))/2
b = m1 + sqrt(12* (m2 - m1.^2))/2

% Exponential
disp("Exponential (lambda):");
lambda = 1 ./ m1

% Erlang
disp("Erlang parameters:");
mean_2 = mean(sortedMatrix .^ 2);
k = round((m1 .^2) ./ (mean_2 - (m1 .^ 2)))
lambda_erlang = k ./ m1
final_lambda_erlang = lambda_erlang(:,1); % perchè per il trace 2 è pari a 0.

% Weibull
% generate randomic values for lambda and k such that the weibull moments
% are similar to the Trace moments
%TRACE 1
m1_perweibull = m1(:,1);
m2_perweibull = m2(:,1);
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibull1 = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibull1);
     % I consider that the max distance between weibull and trace moments is 0.002 
     if abs(wm(1) - m1(:,1)) <= 0.002 & abs(wm(2) - m2(:,1)) <= 0.002
        break;
     end
 end
 disp("Weibull parameters TRACE1:");
 disp(paramweibull1);

% TRACE 2
m1_perweibull = m1(:,2);
m2_perweibull = m2(:,2);
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibull2 = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibull2);
     if abs(wm(1) - m1(:,2)) <= 0.002 & abs(wm(2) - m2(:,2)) <= 0.002
        break;
     end
 end
 disp("Weibull parameters TRACE2:");
 disp(paramweibull2);

 % Pareto
 %TRACE 1
 m1_perpareto = m1(:,1);
 m2_perpareto = m2(:,1);
 while(true)
    pm = [];
    alpha = 0;
    % Because alpha must be > 2
    while alpha <= 2
        alpha = 2 + rand();
    end
    parampareto1 = fsolve(@(x)MM_Pareto_equation(x),[alpha,rand()]);
    pm = pareto_moments(parampareto1);
    if abs(pm(1) - m1(:,1)) <= 0.002 & abs(pm(2) - m2(:,1)) <= 0.002
       break;
    end
 end
 disp("Pareto parameters TRACE1:");
 disp(parampareto1);

 %TRACE 2
 m1_perpareto = m1(:,2);
 m2_perpareto = m2(:,2);
 while(true)
    pm2 = [];
    alpha2 = 0;
    % Because alpha must be > 2
    while alpha2 <= 2
        alpha2 = 2 + rand();
    end
    parampareto2 = fsolve(@(x)MM_Pareto_equation(x),[alpha2,rand()]);
    pm = pareto_moments(parampareto2);
    if abs(pm(1) - m1(:,2)) <= 0.002 & abs(pm(2) - m2(:,2)) <= 0.002
       break;
    end
 end
 disp("Pareto parameters TRACE2:");
 disp(parampareto2);

% Hyper-Exponential
disp("Hyper-Exp parameters:");
column1 = sortedMatrix(:,1);
column2 = sortedMatrix(:,2);
m1_column1 = m1(:,1); 
m1_column2 = m1(:,2); 

%non ha senso calcolare la hyper per il trace 1 perchè cv < 1
paramHyper_t1 = mle(column1, 'pdf', @(column1, l1, l2, p)HyperExp_pdf(column1, [l1, l2, p]), 'start', [0.8/m1_column1, 1.2/m1_column1, 0.4])
paramHyper_t2 = mle(column2, 'pdf', @(column2, l1, l2, p)HyperExp_pdf(column2, [l1, l2, p]), 'start', [0.8/m1_column2, 1.2/m1_column2, 0.4])


% Hypo-Exponential
disp("Hyper-Exp parameters:")
paramHypo_t1 = mle(column1, 'pdf', @(column1, l1, l2)HypoExp_pdf(column1, [l1, l2]), 'start', [1/(0.3*m1_column1), 1/(0.7*m1_column1)])
%non ha senso calcolare la hypo per trace 2 perchè cv > 1
paramHypo_t2 = mle(column2, 'pdf', @(column2, l1, l2)HypoExp_pdf(column2, [l1, l2]), 'start', [1/(0.3*m1_column2), 1/(0.7*m1_column2)])

range = 0:250;
%trace 1
figure;
plot(sortedMatrix(:,1), [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambda(:,1)]), "-", ...
     range, Unif_cdf(range, [a(:,1),b(:,1)]), "-", ...
     range, HypoExp_cdf(range,paramHypo_t1),"-", ...
     range, Erlang_cdf(range, final_lambda_erlang, k(:,1)),"-", ...
     range, Weibull_cdf(range,paramweibull1),"-", ...
     range, Pareto_cdf(range,parampareto1),"-");
title('Trace1');
grid on;
legend('CDF', 'Exp CDF', 'Unif CDF','Hypo CDF','Erlang CDF','Weibull CDF', ...
    'Pareto CDF');

%trace 2
figure;
plot(sortedMatrix(:,2), [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambda(:,2)]), "-", ...
     range, Unif_cdf(range, [a(:,2),b(:,2)]), "-", ...
     range, HyperExp_cdf(range, paramHyper_t2),"-", ...
     range, Weibull_cdf(range,paramweibull2),"-", ...
     range, Pareto_cdf(range,parampareto2),"-");

title('Trace2');
grid on;
legend('CDF', 'Exp CDF', 'Unif CDF','Hyper CDF','Weibull CDF','Pareto CDF');


% FUNCTIONS
function F = weibull_moments(param)
    l = param(1);
    k = param(2);

    F = [];
    F(1) = l * gamma(1 + 1/k);
    F(2) = l^2 * gamma(1 + 2/k);
end

function F = MM_Weibull_equation(param)
	global m1_perweibull;
	global m2_perweibull;
	
	F = weibull_moments(param);
	F(1) = F(1) ./ m1_perweibull - 1;
	F(2) = F(2) ./ m2_perweibull - 1;
end

function F = pareto_moments(p)
	alpha = p(1);
	m = p(2);

	F(1) = (alpha*m)/(alpha-1);
	F(2) = (alpha*(m^2))/(alpha-2);
end

function F = MM_Pareto_equation(p)
	global m1_perpareto;
	global m2_perpareto;

	F = pareto_moments(p);
	F(1) = F(1) ./ m1_perpareto-1;
	F(2) = F(2) ./ m2_perpareto-1;
end

function F = HyperExp_pdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = (x > 0) .* (p1 * l1 * exp(-l1*x) + (1-p1) * l2 * exp(-l2*x));
end

function F = HypoExp_pdf(x, p)
	l1 = p(1);
	l2 = p(2);
	
	F = (x>0) .* (l1*l2/(l1-l2) * (exp(-l2*x) - exp(-l1*x)));
end

function F = Exp_cdf(x, p)
	l = p(1);
	
	F = max(0,1 - exp(-l*x));
end

function F = Unif_cdf(x, p)
	a = p(1);
	b = p(2);
	
	F = max(0, min(1, (x>a) .* (x<b) .* (x - a) / (b - a) + (x >= b)));
end

function F = HypoExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	
	F = (x>0) .* min(1,max(0,1 - l2/(l2-l1) * exp(-l1*x) + l1/(l2-l1) * exp(-l2*x)));
end
function F = HyperExp_cdf(x, p)
	l1 = p(1);
	l2 = p(2);
	p1 = p(3);
	
	F = max(0,1 - p1 * exp(-l1*x) - (1-p1) * exp(-l2*x));
end

function F = Erlang_cdf(x,l,k)
     sumValue = 0;
      for j = 0:(k-1)
        sumValue = sumValue + (l.*x).^j / factorial(j);
      end
      F = 1 - exp(-l.*x) .* sumValue;
    %F = gamcdf(x,k,l);
end

function F = Weibull_cdf(x,p)
	l = p(1);
	k = p(2);

	F = 1 - exp(-(x/l).^k);
end
function F = Pareto_cdf(x,p)
     alpha = p(1);
     m = p(2);
     
     F = [];
     for i=0:(length(x)-1)
         if(i >= m)
            F = [F,1 - (m/i)^alpha];
         else
            F = [F,0];
         end
     end
end














