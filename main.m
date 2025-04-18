clearvars; close all;

numCustomers = 5;
inverters = ones(numCustomers,1);
activeImpedances = [0.00343, 0.04711; 
                    0.00172, 0.02356; 
                    0.00343, 0.04711; 
                    0.00515, 0.07067; 
                    0.00172, 0.02356];

reactiveImpedances = [0.00147, 0.02157;
                      0.00662, 0.09707;
                      0.00147, 0.02157;
                      0.00147, 0.02157;
                      0.00147, 0.02157];

fprintf(repmat('-',1,100));
fprintf('\nSimulating a system with %d customers\n',numCustomers);

sys = PowerSystem(numCustomers,inverters,activeImpedances,reactiveImpedances);

numOutputs = numCustomers; % each customer senses their voltage

