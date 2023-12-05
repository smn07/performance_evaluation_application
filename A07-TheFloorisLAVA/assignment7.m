% Simone Di Ienno, matricola: 225606, codice persona: 10938038 

clear all;
clc;

T_game = 100000;
current_exec = 1;
ns = 0;

dt = 0;
wins = 0;
s = 0;
t = 0;
enter = 1;

while current_exec < T_game
    s = 0;
    enter = 1;
    dt = 0;
    
    while(enter == 1)
	    if s == 0   % i'm in the initial state
            if(rand() < 0.7) %yellow path
                if(rand() < 0.8) %reaching c1
                    ns = 1;
                    dt = -log(prod(rand(4,1))) / 1.5;
                    
                else %player who fall
                    ns = 3; %final point exit
                    dt = -log(rand()) / 0.5;
                end
            else  % light blue path
                if(rand() < 0.3) %reaching c2
                    ns = 2;
                    dt = 3 + (6 - 3) * rand();
                else % players who fall
                    ns = 3; %final point exit
                    dt = -log(rand()) / 0.25;
                end
            end
        end
    
        if s == 1   % i'm in C1
            if(rand() < 0.5) %yellow path
                if(rand() < 0.25) % reaching c2
                    ns = 2;
                    dt = -log(prod(rand(3,1))) / 2;
                else %player who fall
                    ns = 3; %final point exit
                    dt = -log(rand()) / 0.4;
                end
            else  % light path
                if(rand() < 0.6) %reaching c2
                    ns = 2;
                    dt = -log(rand()) / 0.15;
                else % player who fall
                    ns = 3;%final point exit
                    dt = -log(rand()) / 0.2;
                end
            end
        end
    
        if s == 2   % i'm in C2
            if(rand() < 0.6) %green path -> EXIT
                wins = wins + 1;
            end
            dt = -log(prod(rand(5,1))) / 4;
            ns = 3;
        end

        if s == 3 %exit place
            dt = 0;
            ns = 0;
            enter = 0; %condition to exit
        end

        s = ns;
	    t = t + dt;
    end

    current_exec = current_exec + 1;
end

disp("Number of wins:");
disp(wins);
prob_win=wins/T_game;
disp("Probability of winning")
disp(prob_win)
disp("Average duration of a Game")
disp(t/T_game);
disp("Throughput of the system")
finalTime=t+5*T_game;
tHour=finalTime/60;
disp(T_game/tHour);



