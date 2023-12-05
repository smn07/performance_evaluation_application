% Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

m = 2^32;
a = 1664525;
c = 1013904223;
seed = 521191478;
% number of samples to generate
N = 10000;

% for generating values between 0 and 1 -> i'm not considering the seed for
% the first value.
samples = zeros(1, N);
v = seed; 
for i = 1:N
    v = mod(a * v + c, m);
    samples(i) = v / m; 
end

% N1 samples of an exp distrib. with lambda = 0.1
lambdaexp = 0.1;
N1 = 10000;
expsamples = -log(samples(1:N1)) / lambdaexp;

range = 0:25;
% Exp. distribution
sorted_exp_samples = sort(expsamples);
figure;
hold on;
plot(sorted_exp_samples, [1:10000]/10000,'-');
plot(range,Exp_cdf(range,lambdaexp),'-');
title('Exponential');
legend('Generated exponential CDF','Real exponential CDF');
grid on;
hold off;

% cost exp distrib.
cost_exp=calculateCost(10000,expsamples)

% N2 samples of an pareto distrib. with a = 1.5, m = 5
a_pareto = 1.5;
m_pareto = 5;
N2 = 10000;
paretosamples = m_pareto ./ (samples(1:N2)).^(1/a_pareto);

% Pareto distribution
figure;
hold on;
plot(sort(paretosamples), [1:10000]/10000,'-');
plot(range,Pareto_cdf(range,[a_pareto,m_pareto]),'-');
title('Pareto');
legend('Generated Pareto CDF','Real Pareto CDF');
grid on;
hold off;

% cost pareto distrib.
cost_pareto = calculateCost(10000,paretosamples)

% N3 Erlang distribution with k = 4 e lambda = 0.4
k_erlang = 4;
lambda_erlang = 0.4;
N3 = 2500;
X1 = samples(1,1:2500);
X2 = samples(1,2501:5000);
X3 = samples(1,5001:7500);
X4 = samples(1,7501:10000);
X = [X1;X2;X3;X4];
erlangsamples = -log(prod(X)) / lambda_erlang;

% Erlang distribution
figure;
hold on;
plot(sort(erlangsamples), [1:2500]/2500,'-');
plot(range,Erlang_cdf(range,lambda_erlang,k_erlang),'-');
title('Erlang');
legend('Generated Erlang CDF','Real Erlang CDF');
grid on;
hold off;

% cost erlang distrib.
cost_erlang = calculateCost(2500,erlangsamples)

% N4 = 5000 samples of a Hypo-Exponential distribution with rates  lambda1 = 0.5, lambda2 = 0.125
lambda1 = 0.5;
lambda2 = 0.125;
X1 = samples(1,1:5000);
X2 = samples(1,5001:10000);
hypo_expsamples = (-log(X1)/lambda1)-log(X2)/lambda2;

% Hypo-Exp. distribution
figure;
hold on;
plot(sort(hypo_expsamples), [1:5000]/5000,'-');
plot(range,HypoExp_cdf(range,[lambda1,lambda2]),'-');
title('Hypo-Exponential');
legend('Generated Hypo-Exp. CDF','Real Hypo-Exp. CDF');
grid on;
hold off;

% cost hypo-exp distrib.
cost_hypo = calculateCost(5000,hypo_expsamples)

% N5 = 5000 samples of a Hyper-Exponential distribution with rates lambda1 = 0.5, lambda2 = 0.05, p1 = 0.55
l1 = 0.5;
l2 = 0.05;
p1 = 0.55;
X1 = samples(1,1:5000);
X2 = samples(1,5001:10000);

hyper_expsamples = [];

for i = 1:5000
    if X1(i) < p1
        hyper_expsamples(i) = -log(X2(i))/l1;
    else
        hyper_expsamples(i) = -log(X2(i))/l2;
    end
end

% Hyper-Exp. distribution
figure;
hold on;
plot(sort(hyper_expsamples), [1:5000]/5000,'-');
plot(range,HyperExp_cdf(range,[l1,l2,p1]),'-');
title('Hyper-Exponential');
legend('Generated Hyper-Exp. CDF','Real Hyper-Exp. CDF');
grid on;
hold off;

% cost hyper-exp distrib.
cost_hyper = calculateCost(5000,hyper_expsamples)

% FUNCTIONS
function F = Exp_cdf(x, p)
	l = p(1);
	
	F = max(0,1 - exp(-l*x));
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

function F = Erlang_cdf(x,l,k)
     sumValue = 0;
      for j = 0:(k-1)
        sumValue = sumValue + (l.*x).^j / factorial(j);
      end
      F = 1 - exp(-l.*x) .* sumValue;
    %F = gamcdf(x,k,l);
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

function cost = calculateCost(length,vector)
finalSum=0;
 for i=1:length
    if (vector(1,i)) < 10
        finalSum = finalSum + vector(1,i) * 0.01;
    else
        finalSum = finalSum + vector(1,i) * 0.02;
    end
 end
    cost = finalSum;
end






