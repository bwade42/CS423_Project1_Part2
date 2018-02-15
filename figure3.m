function result = figure3()
generations = 10;
N = 10; % total number of elements used in paper 1000 using 10 for testing
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
X = zeros(N,generations);
%%Add Random R values to each element in row
%% This will simulate generation 1 %%
for a = 1: N
   minVal = 3.9;
   maxVal = 4.0;
   r = rand(1) * (maxVal - minVal) + minVal;
   X(a,1) = r;
end

for i = 2: generations
    for j = 1: N
    %%Calculate the Local Dynamics by calling equation #2 %%    
      R = X(j,i-1); %% random fitness value generated earlier %%/
      LocalDynamics = LD(R,i,K);
      
    %%Calculate the Mean Field from equation #4 %%
    %% There is a bug here  you will see when running plot %%
       Mn = IMF(X,N);
       plot(Mn,'*');
       
    %% Plug everything into equation #1 %%
    %% loop through each coupling strength
    %% im just using e5 for now to test the function %%
        result = (1 - e5) * LocalDynamics * (e5 * Mn);
        result = X(j,i); %%update the matrix with the result value %%
        %%plot(result,'*');
    end  
end

end

%% Local Dynamics of each element Equation #2 %%
function localDynamics = LD(R,generation,K)
localDynamics = (R * generation)  * (1 - (generation/K));
end

%% Calculate Global Coupling Strenghth Equation #4 %%
%% X represents the populated matrix %%
function InstaneousMeanField = IMF(X,N)
InstaneousMeanField = sum(X/N);
end