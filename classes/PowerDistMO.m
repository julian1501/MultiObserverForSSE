classdef PowerDistMO
    % PowerDistMO Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numCustomers % Number of customers in the system
        sys % Power Distribution system object
        attack % Attack object
        primaryMO % Large multi observer
        secondaryMO % Small multi observer
    end

    methods
        function obj = PowerDistMO(numCustomers,numAttacks,attackedOutputs,inputFileName)
            % PowerDistMO Construct an instance of this class
            %   Detailed explanation goes here

            fprintf(repmat('-',1,100));
            fprintf('\nCreating a system with %d customers\n',numCustomers);
            
            obj.numCustomers = numCustomers;

            if ~isempty(inputFileName)
                ... % read file extract data
                % Implement when data format is known
            else
                % in case of no file provided use 'standard system'
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
                obj.numCustomers = 5;

                if numAttacks >= obj.numCustomers
                    error("Number of attacks (%d) is larger than the number of customers ()",numAttacks,obj.numCustomers)
                end
            end
            
            obj.sys = PowerSystem(numCustomers,inverters,activeImpedances,reactiveImpedances);
            obj.attack = Attack(obj.numCustomers,numAttacks,attackedOutputs);
            primaryMO = MO(obj.sys,obj.attack)
        
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
        
        
    end
end