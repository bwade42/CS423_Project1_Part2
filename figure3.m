function result = figure3()
generations = 10;
N = 1000; % total number of elements used in paper 1000 using 10 for testing
K = 100; % Carrying capacity

%% Global Coupling Strengths used in Walker Paper %%
e1 = 0.0; %% magenta
e2 = 0.075; %%red
e3 = 0.1; %%blue
e4 = 0.2; %% orange
e5 = 0.25; %% aqua
e6 = 0.3; %% dark green
e7 = 0.3; %% stars
e8 = 0.4; %% black

%% Generate a matrix with elements in rows, and time steps in columns
ModelSystemOutputs = zeros(N,generations);
Y = zeros(N,generations);
%%initialize the system by setting each element at timestep 0 to be 1 %
for a = 1:N
    ModelSystemOutputs(a,1) = 1;
end

for i = 1: generations
    for j = 1: N %% 1 to the total number of elements
        
   
      
      %%Calculate the Local Dynamics by calling equation #2 %% 
      %% this equation depends on equation 1 %%
      LD = LD(i,K);
      
     %%Populate matrix Y with random R variables in rows and, time steps%%
     %% in columns %%
      Y(j,i) = LD;
      
      %%Use equation 4 to determine meanfield for matrix Y %%
      %% This equation depends on equation 2 %%
      MF = IMF(Y,N);
      
      
      %% use equation 1 to calculate output of model system %%
      MS = MS(e1,LD,MF);
        
      
      plot(MS,i ,'*');
       plot(MS2,i ,'+');
     
      
         
    end  
end

end

%% Equation #1 %%
function  ModelSystem = MS(E,LD,MF)
   ModelSystem = ((1 - E) * LD) + (E * MF);
end

%% Local Dynamics of each element Equation #2 %%
function LocalDynamics = LD(generation,K)
    minVal = 3.9;
    maxVal = 4.0;
    r = rand(1) * (maxVal - minVal) + minVal;
LocalDynamics = (r * generation)  * (1 - (generation/K));
end

%% Calculate Global Coupling Strenghth Equation #4 %%
%% X represents the populated matrix %%
function InstaneousMeanField = IMF(X,N)
InstaneousMeanField = sum(X/N);
end



