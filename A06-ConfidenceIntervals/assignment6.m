%Simone Di Ienno, matricola: 225606, codice ID: 10938038

clear all;
clc;

gamma = 0.95; 
alpha = 0.04;
M = 1000; 
d_gamma = norminv((1+gamma)/2);

%SCENARIO 1
%arrival 
l1_hyper = 0.02;
l2_hyper = 0.2;
p1 = 0.1;
Nc=2;

%service times
k_erlang = 10;
l_erlang = 1.5;

% firstly i start with only 1 batch -> and errors are = 1 (value > 0.04)
K1=1;
eU1=1;
eT1=1;
eR1=1;
eJ1=1;
eV1=1;

%vectors of parameters to update once i need a new batch
%mi serve fondamentalmente per calcolare i vari valori che servono per
%calcolare poi gli errori.

U1 = [];%utilization
R1 = [];%response time
J1 = [];%number of jobs
T1 = [];%throughput
V1 = [];%variance of response time

while(eU1 >= alpha || eT1 >= alpha || eR1 >= alpha || eJ1 >= alpha || eV1 >= alpha)
       
    tA1 = 0;
    tC1 = 0;
    Bi1 = 0;
	Wi1 = 0;
	tA01 = 0;
    Bi1 = 0;
	Wi1 = 0;
    vectorR1=[];
    
    % hyper_exp
    samples_A = rand(1,M);
    samples_A2=rand(1,M);

    a1 = [];

    for i = 1:M
        if samples_A(i) < p1
            a1(i) = -log(samples_A2(i))/l1_hyper;
        else
            a1(i) = -log(samples_A2(i))/l2_hyper;
        end
    end
    
    % %Erlang    
    %  samples_S = rand(1,M);
    %  samples_S2=rand(1,M);
    %  samples_S3=rand(1,M);
    %  samples_S4=rand(1,M);
    %  samples_S5=rand(1,M);
    %  samples_S6=rand(1,M);
    %  samples_S7=rand(1,M);
    %  samples_S8=rand(1,M);
    %  samples_S9=rand(1,M);
    %  samples_S10=rand(1,M);
    % 
    %  X=[samples_S;samples_S2;samples_S3;samples_S4;samples_S5;samples_S6;samples_S7;samples_S8;samples_S9;samples_S10];

    %s1 = -log(prod(X)) / l_erlang;


    %ho provato ad usare questa funzione cercando su internet, e dovrebbe
    %venire molto simile rispetto a se li genero con il codice scritto
    %sopra
    s1 = gamrnd(k_erlang, 1/l_erlang, 1000, 1);
    

    for j = 1:M
       if j==1
          tA01=a1(1);
       end  
        tA1 = tA1 + a1(j);
        tC1 = max(tA1, tC1) + s1(j);
		ri1 = tC1 - tA1;
        Bi1 = Bi1 + s1(j);
		Wi1 = Wi1 + ri1;
        vectorR1(j)=ri1;  %devo tener traccia dei singoli ri in modo tale poi da calcolarci la varianza
    end
    
    % values of the current batch with 1000 samples
    Time_1 = tC1 - tA01;
    Ui1 = Bi1 / Time_1;  %utilization
    Ri1 = Wi1 / M;  %response time
    Ji1 = Wi1 / Time_1;  %number of jobs
    Ti1 = Ji1 / Ri1;  %throughput

    % insert new values of this batch into vectors (read above)
    U1(K1) = Ui1;
    R1(K1) = Ri1;
    J1(K1) = Ji1;
    T1(K1) = Ti1;
    V1(K1) = var(vectorR1);

    %now i need to calculate values useful for errors values.
    %utility
    Um1=sum(U1)/K1;
    Us1=sqrt(sum((U1-Um1).^2)/(K1-1));
    if(K1>1) eU1 = 2 * d_gamma * Us1 / sqrt(K1) / Um1; end

    %response
    Rm1=sum(R1)/K1;
    Rs1=sqrt(sum((R1-Rm1).^2)/(K1-1));
    if(K1>1) eR1 = 2 * d_gamma * Rs1 / sqrt(K1) / Rm1; end
    
    %jobs
    Jm1=sum(J1)/K1;
    Js1=sqrt(sum((J1-Jm1).^2)/(K1-1));
    if(K1>1) eJ1 = 2 * d_gamma * Js1 / sqrt(K1) / Jm1; end

    %Throughput
    Tm1=sum(T1)/K1;
    Ts1=sqrt(sum((T1-Tm1).^2)/(K1-1));
    if(K1>1) eT1 = 2 * d_gamma * Ts1 / sqrt(K1) / Tm1; end

    %variance response time
    Vm1=sum(V1)/K1;
    Vs1=sqrt(sum((V1-Vm1).^2)/(K1-1));
    if(K1>1) eV1 = 2 * d_gamma * Vs1 / sqrt(K1) / Vm1; end

    K1=K1+1;
end

% Confidence intervals
cU1=[Um1 - d_gamma * Us1 / sqrt(K1), Um1 + d_gamma * Us1 / sqrt(K1)];
cT1=[Tm1 - d_gamma * Ts1 / sqrt(K1), Tm1 + d_gamma * Ts1 / sqrt(K1)];
cR1=[Rm1 - d_gamma * Rs1 / sqrt(K1), Rm1 + d_gamma * Rs1 / sqrt(K1)];
cJ1=[Jm1 - d_gamma * Js1 / sqrt(K1), Jm1 + d_gamma * Js1 / sqrt(K1)];
cV1=[Vm1 - d_gamma * Vs1 / sqrt(K1), Vm1 + d_gamma * Vs1 / sqrt(K1)];


% Disp results
disp('SCENARIO 1')
disp('Number of iterations needed: ');
disp(K1);
disp('Utilization: ');
disp(cU1);
disp('error for U:')
disp(eU1);
disp('Throughput: ');
disp(cT1);
disp('error for X:')
disp(eT1);
disp('Average number of jobs in the system: ');
disp(cJ1);
disp('error for J:')
disp(eJ1);
disp('Average response time: ');
disp(cR1);
disp('error for R:')
disp(eR1);
disp('Variance of the response time: ');
disp(cV1);
disp('error for V:')
disp(eV1);



%SCENARIO 2
%arrival
l = 0.1;

%service
a = 5;
b = 10;

% firstly i start with only 1 batch -> and errors are = 1 (value > 0.04)
K2=1;
eU2=1;
eT2=1;
eR2=1;
eJ2=1;
eV2=1;

%vectors of parameters to update once i need a new batch
%mi serve fondamentalmente per calcolare i vari valori che servono per
%calcolare poi gli errori.

U2 = [];
R2 = [];
J2 = [];
T2 = [];
V2 = [];

% if all errors are < 0.04 -> accept. 
while(eU2 >= alpha || eT2 >= alpha || eR2 >= alpha || eJ2 >= alpha || eV2 >= alpha)
       
    tA2 = 0;
    tC2 = 0;
    Bi2 = 0;
	Wi2 = 0;
	tA02 = 0;
    Bi2= 0;
	Wi2 = 0;
    vectorR2=[];
    
        %arrival
        samples_A1=rand(1,M);      
        a2 = -log(samples_A1(1:M)) / 0.1;
        %service
        s2 = a + (b - a) * rand(1, M);
    

    for j = 1:M
       if j==1
          tA02=a2(1);
       end  
        tA2 = tA2 + a2(j);
        tC2 = max(tA2, tC2) + s2(j);
        Bi2 = Bi2 + s2(j);
		ri2 = tC2 - tA2;
		Wi2 = Wi2 + ri2;
        vectorR2(j)=ri2;  %devo tener traccia dei singoli ri in modo tale poi da calcolarci la varianza
    end
    
    % values of the current batch with 1000 samples
    Time_2 = tC2 - tA02;
    Ui2 = Bi2 / Time_2;  %utilization
    Ri2 = Wi2 / M;  %response time
    Ji2 = Wi2 / Time_2;  %number of jobs
    Ti2 = Ji2 / Ri2;  %throughput

    % insert new values of this batch into vectors (read above)
    U2(K2) = Ui2;
    R2(K2) = Ri2;
    J2(K2) = Ji2;
    T2(K2) = Ti2;
    V2(K2) = var(vectorR2);

    %now i need to calculate values useful for errors values.
    %utility
    Um2=sum(U2)/K2;
    Us2=sqrt(sum((U2-Um2).^2)/(K2-1));
    if(K2>1) eU2 = 2 * d_gamma * Us2 / sqrt(K2) / Um2; end

    %response
    Rm2=sum(R2)/K2;
    Rs2=sqrt(sum((R2-Rm2).^2)/(K2-1));
    if(K2>1) eR2 = 2 * d_gamma * Rs2 / sqrt(K2) / Rm2; end
    
    %jobs
    Jm2=sum(J2)/K2;
    Js2=sqrt(sum((J2-Jm2).^2)/(K2-1));
    if(K2>1) eJ2 = 2 * d_gamma * Js2 / sqrt(K2) / Jm2; end

    %Throughput
    Tm2=sum(T2)/K2;
    Ts2=sqrt(sum((T2-Tm2).^2)/(K2-1));
    if(K2>1) eT2 = 2 * d_gamma * Ts2 / sqrt(K2) / Tm2; end

    %variance response time
    Vm2=sum(V2)/K2;
    Vs2=sqrt(sum((V2-Vm2).^2)/(K2-1));
    if(K2>1) eV2 = 2 * d_gamma * Vs2 / sqrt(K2) / Vm2; end

    K2=K2+1;
end

% Confidence intervals
cU2=[Um2 - d_gamma * Us2 / sqrt(K2), Um2 + d_gamma * Us2 / sqrt(K2)];
cT2=[Tm2 - d_gamma * Ts2 / sqrt(K2), Tm2 + d_gamma * Ts2 / sqrt(K2)];
cR2=[Rm2 - d_gamma * Rs2 / sqrt(K2), Rm2 + d_gamma * Rs2 / sqrt(K2)];
cJ2=[Jm2 - d_gamma * Js2 / sqrt(K2), Jm2 + d_gamma * Js2 / sqrt(K2)];
cV2=[Vm2 - d_gamma * Vs2 / sqrt(K2), Vm2 + d_gamma * Vs2 / sqrt(K2)];


% Disp results
disp('SCENARIO 2')
disp('Number of iterations needed: ');
disp(K2);
disp('Utilization: ');
disp(cU2);
disp('error for U:')
disp(eU2);
disp('Throughput: ');
disp(cT2);
disp('error for X:')
disp(eT2);
disp('Average number of jobs in the system: ');
disp(cJ2);
disp('error for J:')
disp(eJ2);
disp('Average response time: ');
disp(cR2);
disp('error for R:')
disp(eR2);
disp('Variance of the response time: ');
disp(cV2);
disp('error for V:')
disp(eV2);