%% Simone Di Ienno, matricola: 225606, codice ID: 10938038
% PROJECT D
clear all;
clc;

%% Global variables for weibull and pareto equation functions
global m1_perweibull;
global m2_perweibull;
global m1_perpareto;
global m2_perpareto;

% Times are expressed in days
% I know the distribution for writing and shooting stages (Erlang)

%% Fitting for Audio Editing stage
audioTrace = importdata('TraceD-A.txt');
sortedAudioTrace = sort(audioTrace);
N = length(audioTrace);
range = [0:0.1:max(audioTrace)];

moment1 = mean(audioTrace);
moment2 = sum(audioTrace .^ 2) / N;

%Exponential
lambdaExpAudio = 1 / moment1;

%Uniform
a = moment1 - sqrt(12* (moment2 - moment1 ^ 2)) / 2;
b = moment1 + sqrt(12* (moment2 - moment1 ^ 2)) / 2;

%Weibull
m1_perweibull = moment1;
m2_perweibull = moment2;
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibullAudio = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibullAudio);
     % I consider that the max distance between weibull and trace moments is 0.002 
     if abs(wm(1) - moment1) <= 0.002 & abs(wm(2) - moment2) <= 0.002
        break;
     end
 end

 %Pareto
 m1_perpareto = moment1;
 m2_perpareto = moment2;
 pm = [];
 while(true)
    pm = [];
    alpha = 0;
    % Because alpha must be > 2
    while alpha <= 2
        alpha = 2 + rand();
    end
    paramparetoAudio = fsolve(@(x)MM_Pareto_equation(x),[alpha,rand()]);
    pm = pareto_moments(paramparetoAudio);
    if abs(pm(1) - moment1) <= 0.002 & abs(pm(2) - moment2) <= 0.002
       break;
    end
 end

%Hyper-Exponential
paramHyperAudio = mle(sortedAudioTrace, 'pdf', @(sortedAudioTrace, l1, l2, p)HyperExp_pdf(sortedAudioTrace, [l1, l2, p]), 'start', [0.8/moment1, 1.2/moment1, 0.4]);

% PRINT VALUES
disp('PARAMETRI HYPEREXP AUDIO EDITING:');
disp(paramHyperAudio);
disp('PARAMETRI EXP AUDIO EDITING:');
disp(lambdaExpAudio);
disp('PARAMETRI WEIBULL AUDIO EDITING:');
disp(paramweibullAudio);
disp('PARAMETRI PARETO AUDIO EDITING:');
disp(paramparetoAudio);
disp('PARAMETRI UNIFORM AUDIO EDITING:');
disp([a,b]);

figure;
plot(sort(audioTrace), [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambdaExpAudio]), "-", ...
     range, Unif_cdf(range, [a,b]), "-", ...
     range, Weibull_cdf(range,paramweibullAudio),"-", ...
     range, Pareto_cdf(range,paramparetoAudio),"-", ...
     range, HyperExp_cdf(range, paramHyperAudio),"-");
title('Audio Editing trace');
xlim([0,50]);  % Set the x-axis limits
grid on;
legend('AudioTrace', 'Exp CDF', 'Unif CDF','Weibull CDF', ...
    'Pareto CDF', 'Hyper CDF');

%% Fitting for Video Editing stage
videoTrace = importdata('TraceD-V.txt');
sortedVideoTrace1 = sort(videoTrace);
sortedVideoTrace = sortedVideoTrace1(2:end);
N = length(sortedVideoTrace);
range = [0:0.1:max(videoTrace)];

moment1 = mean(videoTrace);
moment2 = sum(videoTrace .^ 2) / N;

%Exponential
lambdaExpVideo = 1 / moment1;

%Uniform
a = moment1 - sqrt(12* (moment2 - moment1 ^ 2)) / 2;
b = moment1 + sqrt(12* (moment2 - moment1 ^ 2)) / 2;

%Weibull
m1_perweibull = moment1;
m2_perweibull = moment2;
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibullVideo = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibullVideo);
     % I consider that the max distance between weibull and trace moments is 0.002 
     if abs(wm(1) - moment1) <= 0.002 & abs(wm(2) - moment2) <= 0.002
        break;
     end
 end

 %Pareto
 m1_perpareto = moment1;
 m2_perpareto = moment2;
 pm = [];
 while(true)
    pm = [];
    alpha = 0;
    % Because alpha must be > 2
    while alpha <= 2
        alpha = 2 + rand();
    end
    paramparetoVideo = fsolve(@(x)MM_Pareto_equation(x),[alpha,rand()]);
    pm = pareto_moments(paramparetoVideo);
    if abs(pm(1) - moment1) <= 0.002 & abs(pm(2) - moment2) <= 0.002
       break;
    end
 end

%Hyper-Exponential
paramHyperVideo = mle(sortedVideoTrace, 'pdf', @(sortedVideoTrace, l1, l2, p)HyperExp_pdf(sortedVideoTrace, [l1, l2, p]), 'start', [0.8/moment1, 1.2/moment1, 0.4]);

% PRINT VALUES
disp('PARAMETRI HYPEREXP VIDEO EDITING:');
disp(paramHyperVideo);
disp('PARAMETRI EXP VIDEO EDITING:');
disp(lambdaExpVideo);
disp('PARAMETRI WEIBULL VIDEO EDITING:');
disp(paramweibullVideo);
disp('PARAMETRI PARETO VIDEO EDITING:');
disp(paramparetoVideo);
disp('PARAMETRI UNIFORM VIDEO EDITING:');
disp([a,b]);

figure;
plot(sortedVideoTrace, [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambdaExpVideo]), "-", ...
     range, Unif_cdf(range, [a,b]), "-", ...
     range, Weibull_cdf(range,paramweibullVideo),"-", ...
     range, Pareto_cdf(range,paramparetoVideo),"-", ...
     range, HyperExp_cdf(range, paramHyperVideo),"-");
title('Video Editing trace');
xlim([0,50]);  % Set the x-axis limits
grid on;
legend('VideoTrace', 'Exp CDF', 'Unif CDF','Weibull CDF', ...
    'Pareto CDF', 'Hyper CDF');

%% Fitting VFX
compositingTrace = importdata('TraceD-X.txt');
sortedVfxTrace = sort(compositingTrace);
N = length(sortedVfxTrace);
range = [0:0.1:max(compositingTrace)];

moment1 = mean(compositingTrace);
moment2 = sum(compositingTrace .^ 2) / N;

%Exponential
lambdaExpVfx = 1 / moment1;

%Uniform
a = moment1 - sqrt(12* (moment2 - moment1 ^ 2)) / 2;
b = moment1 + sqrt(12* (moment2 - moment1 ^ 2)) / 2;

%Weibull
m1_perweibull = moment1;
m2_perweibull = moment2;
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibullVfx = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibullVfx);
     % I consider that the max distance between weibull and trace moments is 0.002 
     if abs(wm(1) - moment1) <= 0.002 & abs(wm(2) - moment2) <= 0.002
        break;
     end
 end

 %Pareto
 m1_perpareto = moment1;
 m2_perpareto = moment2;
 pm = [];
 while(true)
    pm = [];
    alpha = 0;
    % Because alpha must be > 2
    while alpha <= 2
        alpha = 2 + rand();
    end
    paramparetoVfx = fsolve(@(x)MM_Pareto_equation(x),[alpha,rand()]);
    pm = pareto_moments(paramparetoVfx);
    if abs(pm(1) - moment1) <= 0.002 & abs(pm(2) - moment2) <= 0.002
       break;
    end
 end

%Hyper-Exponential
paramHyperVfx = mle(sortedVfxTrace, 'pdf', @(sortedVfxTrace, l1, l2, p)HyperExp_pdf(sortedVfxTrace, [l1, l2, p]), 'start', [0.8/moment1, 1.2/moment1, 0.4]);

% PRINT VALUES
disp('PARAMETRI HYPEREXP VFX:');
disp(paramHyperVfx);
disp('PARAMETRI EXP VFX:');
disp(lambdaExpVfx);
disp('PARAMETRI WEIBULL VFX:');
disp(paramweibullVfx);
disp('PARAMETRI PARETO VFX:');
disp(paramparetoVfx);
disp('PARAMETRI UNIFORM VFX:');
disp([a,b]);

figure;
plot(sortedVfxTrace, [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambdaExpVfx]), "-", ...
     range, Unif_cdf(range, [a,b]), "-", ...
     range, Weibull_cdf(range,paramweibullVfx),"-", ...
     range, Pareto_cdf(range,paramparetoVfx),"-", ...
     range, HyperExp_cdf(range, paramHyperVfx),"-");
title('VFX trace');
xlim([0,50]);  % Set the x-axis limits
grid on;
legend('VfxTrace', 'Exp CDF', 'Unif CDF','Weibull CDF', ...
    'Pareto CDF', 'Hyper CDF');

%% Fitting Compositing
compositingTrace = importdata('TraceD-C.txt');
sortedCompositingTrace = sort(compositingTrace);
N = length(sortedCompositingTrace);
range = [0:0.1:max(compositingTrace)];

moment1 = mean(compositingTrace);
moment2 = sum(compositingTrace .^ 2) / N;

%Exponential
lambdaExpCompositing = 1 / moment1;

%Uniform
a = moment1 - sqrt(12* (moment2 - moment1 ^ 2)) / 2;
b = moment1 + sqrt(12* (moment2 - moment1 ^ 2)) / 2;

%Weibull
m1_perweibull = moment1;
m2_perweibull = moment2;
wm = [];
 while(true)
     wm = [];
     ar = rand();
     br = rand();
     paramweibullCompositing = fsolve(@(x)MM_Weibull_equation(x),[ar,br]);
     wm = weibull_moments(paramweibullCompositing);
     % I consider that the max distance between weibull and trace moments is 0.002 
     if abs(wm(1) - moment1) <= 0.002 & abs(wm(2) - moment2) <= 0.002
        break;
     end
 end

 %Pareto
 m1_perpareto = moment1;
 m2_perpareto = moment2;
 pm = [];
 while(true)
    pm = [];
    alpha = 0;
    % Because alpha must be > 2
    while alpha <= 2
        alpha = 2 + rand();
    end
    paramparetoCompositing = fsolve(@(x)MM_Pareto_equation(x),[alpha,rand()]);
    pm = pareto_moments(paramparetoCompositing);
    if abs(pm(1) - moment1) <= 0.002 & abs(pm(2) - moment2) <= 0.002
       break;
    end
 end

%Hyper-Exponential
paramHyperCompositing = mle(sortedCompositingTrace, 'pdf', @(sortedCompositingTrace, l1, l2, p)HyperExp_pdf(sortedCompositingTrace, [l1, l2, p]), 'start', [0.8/moment1, 1.2/moment1, 0.4], 'Options', statset('MaxFunEvals', 1e5, 'MaxIter', 1e5));

% PRINT VALUES
disp('PARAMETRI HYPEREXP COMPOSITING:');
disp(paramHyperCompositing);
disp('PARAMETRI EXP COMPOSITING:');
disp(lambdaExpCompositing);
disp('PARAMETRI WEIBULL COMPOSITING:');
disp(paramweibullCompositing);
disp('PARAMETRI PARETO COMPOSITING:');
disp(paramparetoCompositing);
disp('PARAMETRI UNIFORM COMPOSITING:');
disp([a,b]);

figure;
plot(sortedCompositingTrace, [1:N]/N, ".", ...
     range, Exp_cdf(range, [lambdaExpCompositing]), "-", ...
     range, Unif_cdf(range, [a,b]), "-", ...
     range, Weibull_cdf(range,paramweibullCompositing),"-", ...
     range, Pareto_cdf(range,paramparetoCompositing),"-", ...
     range, HyperExp_cdf(range, paramHyperCompositing),"-");
title('Compositing trace');
xlim([0,50]);  % Set the x-axis limits
grid on;
legend('CompositingTrace', 'Exp CDF', 'Unif CDF','Weibull CDF', ...
    'Pareto CDF', 'Hyper CDF');

%% FUNCTIONS
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

function F = Exp_cdf(x, p)
	l = p(1);
	
	F = max(0,1 - exp(-l*x));
end

function F = Unif_cdf(x, p)
	a = p(1);
	b = p(2);
	
	F = max(0, min(1, (x>a) .* (x<b) .* (x - a) / (b - a) + (x >= b)));
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
end

function F = Weibull_cdf(x,p)
	l = p(1);
	k = p(2);

	F = 1 - exp(-(x/l).^k);
end

function F = Pareto_cdf(x,p)
     alpha = p(1);
     m = p(2);
     
     F = (1- (m ./ x) .^ alpha) .* (x >= m);
end

