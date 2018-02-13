%% E = global coupling strength return result from mn
%% N = total number of elements
%% i = element
%% n = current timestep(generation)
%% R = Reprodutive fitness of population i
%% K = carrying capacity

function result = figure3()

R = 3.9; % Should be a random value between [3.9,4.0]
N = 100; % total number of elements
K = 100; % Carrying capacity
E = 1.0;
element1 = 1;
genrations = 1000;

hold on
for i = 1: genrations
LocalDynamics = LD(R,element1,i,K); %% fi(xi*n)

%%Mn = MF(LocalDynamics, N); %% Equation #3
mn = IMF(LocalDynamics,N); %% Equation #4

result = (1 - E) * LocalDynamics + (mn * E);

plot(i,result,'*')
end
hold off
end

function LocalDynamics = LD(R,element,generation,K)
LocalDynamics = (R * element * generation)  * (1 - (element * generation) / K);
end

function MeanField = MF(LD, N)
MeanField = LD * N;
end
%% Calculate Global Coupling Strenghth %%
function InstaneousMeanField = IMF(LD,N)
syms k
InstaneousMeanField = symsum(LD * k, k, 1, N);
end
