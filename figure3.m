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
XtoMtransferEntropy = zeros(50,1);
MtoXtransferEntropy = zeros(50,1);
counter = 1;

for e = 0:0.02:1
    
    nextCouplingValue = e

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
    
    
    randomlySelectedPopulation = round((rand(1)*N))';
    XforTE = zeros(generations,1);
    MforTE = zeros(generations,1);
    genCounter = 1;
    
    for t=1:generations
        nextX = X(randomlySelectedPopulation,t);
        nextM = instMeanField(t);
        nextDiscreteX = int64(round(nextX*100)');
        nextDiscreteM = int64(round(nextM*100)');
        XforTE(genCounter) = nextDiscreteX;
        MforTE(genCounter) = nextDiscreteM;
        genCounter = genCounter +1;
    end
    
    %disp(XforTE);
    %disp(MforTE);
    
    %Calculate transfer entropy
    XtoMtransferEntropy(counter) = TE(XforTE, MforTE);
    MtoXtransferEntropy(counter) = TE(MforTE, XforTE);
    
    counter = counter + 1; 
    
end    
    
    eValues = (0:0.02:1);
    plot(eValues, XtoMtransferEntropy,[0 0 0],'--',eValues,MtoXtransferEntropy,[0 0 0]);

% Local Dynamics of each element Equation #2 %%
    function localDynamics = LD(rValue,population)
        localDynamics = (rValue * population)  * (1 - (population/K));

        if localDynamics < 0
            localDynamics = 0;
        end

    end

    function transferEntropy = TE(sourceData, destinationData)
        
        % Add JIDT jar library to the path, and disable warnings that it's already there:
        warning('off','MATLAB:Java:DuplicateClass');
        javaaddpath('C:\Users\johna\Desktop\JIDT\infodynamics.jar');
        % Add utilities to the path
        addpath('C:\Users\johna\Desktop\JIDT\demos\octave');

        % 0. Load/prepare the data:
        %data = load('C:\Users\johna\Documents\MATLAB\nextXtoM.txt');
        % Column indices start from 1 in Matlab:
        source = octaveToJavaDoubleArray(sourceData);
        destination = octaveToJavaDoubleArray(destinationData);

        % 1. Construct the calculator:
        calc = javaObject('infodynamics.measures.continuous.kernel.TransferEntropyCalculatorKernel');
        % 2. Set any properties to non-default values:
        % No properties were set to non-default values
        % 3. Initialise the calculator for (re-)use:
        calc.initialise();
        % 4. Supply the sample data:
        calc.setObservations(source, destination);
        % 5. Compute the estimate:
        result = calc.computeAverageLocalOfObservations();

        fprintf('TE_Kernel(col_0 -> col_1) = %.4f bits\n', ...
            result);

        transferEntropy = result
               
    end

end
