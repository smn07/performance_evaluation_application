% SIMONE DI IENNO, MATRICOLA: 225606, CODICE ID: 10938038

clear all;

fileID = fopen('barrier.log', 'r');

% Arrays IN and OUT
IN = [];
OUT = [];

% read the first line of the barrier.log file
currentLine = fgetl(fileID);
while ischar(currentLine)
    % value between the first []
    tokens = regexp(currentLine, '\[(.*?)\]', 'tokens');
    if ~isempty(tokens)
        stringa = tokens{1}{1};

        % Extract the tenth of second part for organizing the datetime also
        % with tenth of seconds
        tenthOfSecond = str2double(extractAfter(stringa, 18));
        dt = datetime(stringa(1:18), 'InputFormat', 'uuuu:DDD:HH:mm:ss:');
        dt = dt + seconds(tenthOfSecond * 0.1); % i'm adding tenthOfSecond Value

        dt.Format = 'dd-MMM-yyyy HH:mm:ss.S'; % format for displaying also tenthOfSecondValue
        
        % Distinguish between IN and OUT
        if contains(currentLine, 'IN')
            IN = [IN; dt];
        else
            OUT = [OUT; dt];
        end
    end
    
    % Read the next line
    currentLine = fgetl(fileID);
end

fclose(fileID);

%let's start to calculate T
diff = OUT(end)-IN(1);

T = seconds(diff);

A = length(IN); % #arrivals
C = length(OUT);% #completions

% Arrival rate and throughput
ArrivalRate = A/T;
Throughput = C/T;

% Average inter-arrival time
a = IN(2:end) - IN(1:end-1);
averageInterArrivalTime = seconds(sum(a)/A);
% Or i could calculate it with this formula:
averageInterArrivalTime = 1/ArrivalRate;

% Service time i-th
si = OUT(1) - IN(1);
for i = 2:length(OUT)
    si = [si;OUT(i) - max(OUT(i-1), IN(i))];
end

% Average Service time
S = seconds(sum(si)/length(OUT));
% Utilization
U = Throughput*S;

% Average number of jobs
ri = OUT - IN; % i-th response time
W = sum(ri);
ANJ = seconds(W/T); %average number of jobs

% Average Response time
ART = seconds(W/length(OUT));

%Pr(N = 0)
INv = [seconds(IN-datetime('0000-01-01','InputFormat','uuuu-MM-dd')), ones(length(IN),1)];
OUTv = [seconds(OUT-datetime('0000-01-01','InputFormat','uuuu-MM-dd')), -ones(length(OUT),1)];

matrix = [INv; OUTv]; % matrix with ones and -ones

matrix = sortrows(matrix,1); % sort with respect to the first column
matrix(:,3) = cumsum(matrix(:,2)); % cumulative sum
matrix(1:end-1, 4) = matrix(2:end,1) - matrix(1:end-1,1);

% Prob of having N=0
matrix(:,5) = (matrix(:,3) == 0) .* matrix(:,4);
PN0 = sum(matrix(:, 5))/T;

% Prob of having N=1
matrix(:,6) = (matrix(:,3) == 1) .* matrix(:,4);
PN1 = sum(matrix(:, 6))/T;

% Prob of having N=2
matrix(:,7) = (matrix(:,3) == 2) .* matrix(:,4);
PN2 = sum(matrix(:, 7))/T;

%Prob response time < 30 secs
PR30 = sum(seconds(ri) < 30)/length(OUT);

%Prob response time < 3 min
PR180 = sum(seconds(ri) < 180)/length(OUT);

%Prob inter-arrival time < 1 min
PI60 = sum(seconds(a) < 60)/length(OUT);

%Probab service time larger than 1 min
PS60 = sum(seconds(si) > 60)/length(OUT);

%PRINT VALUES FOR THE ASSIGNMENT
disp('Arrival Rate =')
disp(ArrivalRate)
disp('Throughput =')
disp(Throughput)
disp('Average inter-arrival time =')
disp(averageInterArrivalTime)
disp('Utilization =')
disp(U)
disp('Average Service Time =')
disp(S)
disp('Average Number of Jobs =')
disp(ANJ)
disp('Average Response Time =')
disp(ART)
disp('Probability 0 parts in the machine =')
disp(PN0)
disp('Probability 1 parts in the machine =')
disp(PN1)
disp('Probability 2 parts in the machine =')
disp(PN2)
disp('Probability response time < 30 secs =')
disp(PR30)
disp('Probability response time < 3 min =')
disp(PR180)
disp('Probability inter-arrival time < 1 min =')
disp(PI60)
disp('Probability servie time > 1 min =')
disp(PS60)

















