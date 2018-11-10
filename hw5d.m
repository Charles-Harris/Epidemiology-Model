function hw5d

clc;

% Number of martians
N = 1000;

% Disease incubation variables
a1 = 4;
a2 = 9;

% Recovery Variables
b1 = 13;
b2 = 15;

% Immunity Variables
c1 = 6;
c2 = 16;

% Probability of Disease Transfer
DT = 0.6;

% Number of Healthy Martians
NH = 999;

% Number of Infected Martians
NI = 0;

% Number of Contaminated Martians
NC = 1;

% Number of Healthy Martians with Immunity
NHI = 0;

% Number of Days
ND = 500;

% Number of Simulations
NSIM = 100;

% 3D array for storage (Number of Days by number of possible outcomes by number of
% simulations)

Notepad = zeros(ND,4,NSIM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation Loop
for simulationNumber = 1:NSIM
    
    %1. Initialize Subpopulations
    
    % Healthy and Not Immune
    m1 = zeros(NH,1);
    im1 = zeros(NH,1);

    % Infected
    m2 = -floor((a2-a1+1)*rand(NI,1))-a1;
    im2 = zeros(NI,1);

    % Contagious
    m3 =  floor((b2-b1+1)*rand(NC,1))+b1;
    im3 = zeros(NC,1);

    % Healthy and Immune
    m4 = zeros(NHI, 1);
    im4 = -floor((c2-c1+1)*rand(NHI,1))-c1;

    m = [m1' m2' m3' m4']';
    im = [im1' im2' im3' im4']';
    
    Notepad(1,:,simulationNumber) = [NI, NC, NHI, NH];


    %2. Sample Period Loop

    % Exposure Loop Begins
    for day =2:ND
        % a. Randomly Expose population
        for infection = 1:(N/2)
            dt = floor(rand+DT);
            if dt == 1;
                % Algorithm for encounters between martians
                I = floor(N*rand)+1;
                J = floor(N*rand)+1;
                % Martian I and Martian J cant be the same Martian
                while I == J
                    J = floor(N*rand) + 1;
                end
                % Decides whether to infect Martian I or Martian J
                if m(I) > 0 && m(J) == 0 && im(J) >=0;
                   m(J) = -floor((a2-a1+1)*rand)-a1;
                end
                if m(J) > 0 && m(I) == 0 && im(I) >= 0;
                    m(I) = - floor((a2-a1+1)*rand) - a1;
                end
            end
        end
        
        % b. Record NI, NC, NHI, NH
        
        % Infected
        Notepad(day,1,simulationNumber) = sum(m<0);
          

        % Contagious
        Notepad(day,2, simulationNumber) = sum(m>0);
          
        
        % Healthy and Immune
        Notepad(day,3, simulationNumber) = sum (m==0 & im < 0);
         

        % Healthy and Not Immune
        Notepad(day,4, simulationNumber) = sum(m==0 & im >=0);
   
        
         % c. Update Health and Disease Status
          for index = 1:N
          if m(index) == 0
              im(index) = im(index) + 1;
          end
          
          if m(index) > 0
             m(index) = m(index)-1;
             if m(index) == 0
                 im(index) = - floor((c2-c1+1)*rand)-c1;
             end
          end
          
          if m(index) < 0
              m(index)= m(index) + 1;
              if m(index)== 0
                  m(index) = floor((b2-b1+1)*rand)+b1;
              end
          end
          end
    end
end

% Plotting

% Plot One
figure;
hold on;
for iteration = 1:NSIM
    hold on
    %subplot(2,2,1);
    plotOne = plot(Notepad(:,1,iteration),'y');
end
title('Infected People Superposition');
ylabel('Number of People');
xlabel('Number of Days Passed');
legend(plotOne,'NI Population');
axis([1 ND 0 N]);

% Plot Two
figure;
hold on;
for iteration = 1:NSIM
    hold on
    %subplot(2,2,2);
    plotTwo = plot(Notepad(:,2,iteration),'r');
end
title('Contagious People Superposition');
ylabel('Number of People');
xlabel('Number of Days Passed');
legend(plotTwo,'NC Population');
axis([1 ND 0 N]);

%Plot 3
figure;
hold on;
for iteration = 1:NSIM
    hold on
    %subplot(2,2,3);
    plotThree = plot(Notepad(:,3,iteration),'k');
end
title('Helathy People with Immunity Superposition');
ylabel('Number of People');
xlabel('Number of Days Passed');
legend(plotThree,'NHI Population');
axis([1 ND 0 N]);

% Plot 4
figure;
hold on;
for iteration = 1:NSIM
    hold on
    % subplot(2,2,4);
    plotFour = plot(Notepad(:,4,iteration),'G');
end
title('Healthy People without Immunity Superposition');
ylabel('Number of People');
xlabel('Number of Days Passed');
legend(plotFour,'NH Population');
axis([1 ND 0 N]);

%Plotting the histograms

% Summing totals for each simulation
for simulation = 1:NSIM
    Infection(simulation) = sum(Notepad(ND,1,simulation));
    Contagion(simulation) = sum(Notepad(ND,2,simulation));
    Immunity(simulation) = sum(Notepad(ND,3,simulation));
    Healthy(simulation) = sum(Notepad(ND,4,simulation));
end
    



% Establishing BinWidth

binWidthVector = 0:10:1000;

%Histogram 1
%subplot(2,2,1); 
figure;
hist(Infection,binWidthVector);
ylabel('Num Sims');
xlabel('Num People');
title('Infectious People');
xlim([0 N]);


%Histogram 2
%subplot(2,2,2);
figure; 
hist(Contagion,binWidthVector);
ylabel('Num Sims');
xlabel('Num People');
title('Contagious People');
xlim([0 N]);



%Histogram 3
%subplot(2,2,3); 
figure;
hist(Immunity,binWidthVector);
ylabel('Num Sims');
xlabel('Num People');
title('Immune People');
xlim([0 N]);



%Histogram 4
%subplot(2,2,4);
figure;
hist(Healthy,binWidthVector);
ylabel('Num Sims');
xlabel('Num People');
title('Healthy People (No immunity)');
xlim([0 N]);
end