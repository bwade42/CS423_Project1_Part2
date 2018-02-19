function result = figure3()
generations = 10000;
N = 1000; % total number of elements used in paper 1000 using 10 for testing
K = 100; % Carrying capacity

%% Global Coupling Strengths used in Walker Paper %%
e1 = 0.0; %% magenta
e2 = 0.075; %%red
e3 = 0.1; %%blue
e4 = 0.2; %% orange
e5 = 0.225; %% aqua
e6 = 0.25; %% dark green
e7 = 0.3; %% stars
e8 = 0.4; %% black

globalCouplingValues = [e1 e2 e3 e4 e5 e6 e7 e8];

for e=1:8
    
    nextCouplingValue = globalCouplingValues(e);

    %% Generate a matrix with elements in rows, and time steps in columns
    X = zeros(N,generations);
    rValues = zeros(N);
    localDynamics = zeros(N,generations);
    instMeanField = zeros(generations);
    instDyanmicsOfMeanField = zeros(generations);

    %%Generate a random fitness parameter 'r' for each population
    for a = 1: N
       minVal = 3.9;
       maxVal = 4.0;
       r = rand(1) * (maxVal - minVal) + minVal;
       rValues(a) = r;

       %%Initialize x0 = 1
       X(a,1) = 1;
       localDynamics(a,1) = 0;
       instMeanField(1) = 1;
       instDyanmicsOfMeanField(1) = 1;
    end

    for i = 2: generations
        %%Calculate the local dynamics for this generation
        for j = 1: N
        %%Calculate the Local Dynamics by calling equation #2 %%    
          localDynamics(j,i) = LD(rValues(j),X(j,i-1));
        end

        %%Use the local dynamics of this generation to calculate the
        %%Instantaneous dynamics of the mean-field
        for j = 1: N
            %%Add the local dynamics to the instantaneous dynamics of the mean field for this generation%%
            instDyanmicsOfMeanField(i) = instDyanmicsOfMeanField(i) + localDynamics(j,i);
        end

        instDyanmicsOfMeanField(i) = instDyanmicsOfMeanField(i)/N;

        %%Finally, calculate the x value for this generation
        for j = 1: N
            X(j,i) = (1-nextCouplingValue)*localDynamics(j,i) + (nextCouplingValue*instDyanmicsOfMeanField(i));
        end
    end

    %%Now calculate the value of the instantaneous mean field for each
    %%generation
    for i= 2: generations
        for j= 1: N
           instMeanField(i) = instMeanField(i) + X(j,i);
        end

        instMeanField(i) = instMeanField(i)/N;
    end

    meanFieldX = instMeanField(1:N-1);
    meanFieldY = instMeanField(2:N);
    
    switch e
        case 1
            nextGraphColor = [1 0 1];
        case 2
            nextGraphColor = [1 0 0];
        case 3
            nextGraphColor = [0 0 1];
        case 4
            nextGraphColor = [0.9 0.6 0.1];
        case 5
            nextGraphColor = [0 1 1];
        case 6
            nextGraphColor = [0 0.39 0];
        case 7
            nextGraphColor = [0.2 0 0.4];
        case 8
            nextGraphColor = [0 0 0];       
    end

    scatter(meanFieldX, meanFieldY, 4,nextGraphColor,'.');
    hold on


    %
end

% Local Dynamics of each element Equation #2 %%
    function localDynamics = LD(rValue,population)
        localDynamics = (rValue * population)  * (1 - (population/K));

        if localDynamics < 0
            localDynamics = 0;
        end

    end

end
