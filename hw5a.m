function hw5a

clc;

% Number of martians
N = 100;

% Disease incubation variables
a1 = 5;
a2 = 10;

% Recovery Variables
b1 = 10000;
b2 = 20000;

% Immunity Variables
c1 = 0;
c2 = 0;

% Probability of Disease Transfer
DT = 0.6;

% Number of Healthy Martians
NH = 99;

% Number of Infected Martians
NI = 0;

% Number of Contaminated Martians
NC = 1;

% Number of Healthy Martians with Immunity
NHI = 0;

% Number of Days
ND = 60;

% Number of Simulations
NSIM = 5;

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
figure;
hold on;
for simu = 1:NSIM
    plot(Notepad(:,1,simu),'y');
    plot(Notepad(:,2,simu),'r');
    plot(Notepad(:,4,simu),'k');
end
ylabel('Number of Martians');
xlabel('Number of Days Passed');
title('Number of Healthy w/0 Immunity (NH), Infected(NI), and Contagious (NC)');
legend('Infected(NI)', 'Contagious (NC)','Healthy w/o Immunity (NH)');
axis ([1 ND 0 N]);

% Nakes the table
X = ceil(rand()*NSIM);
A = [(1:1:ND)' Notepad(:,3,X) Notepad(:,4,X) Notepad(:,1,X) Notepad(:,2,X)];
T = table(A(:,1), A(:,2), A(:,3), A(:,4), A(:,5), 'VariableNames', {'ND' 'NHI' 'NH' 'NI' 'NC'});
disp(T);
end




        
            
        
        
    
    
    

    
        
        
    