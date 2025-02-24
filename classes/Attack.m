classdef Attack
    % Attack Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numAttacks
        attackList
    end

    methods
        function obj = Attack(numOutputs,numAttacks,attackedOutputs)
            % Attack Construct an instance of this class
            %   Detailed explanation goes here
            obj.numAttacks = numAttacks;
            if ~ (numOutputs > 2* numAttacks)
               warning('The number of outputs is not larger then twice the number of attacked outputs %3.0f <= %3.0f',numOutputs,numAttacks); 
            end

            % pick random outputs to attack if none have been specified
            if size(attackedOutputs) == [0,0]
                attackedOutputs = sort(randperm(numCustomers,numAttacks));
            end
            fprintf("The attacked outputs are: \n");
            disp(attackedOutputs);
            
            % loop over outputs and set booleans indicating if it is
            % attacked
            attackList = zeros(numOutputs,1);
            for i = 1:1:obj.numAttacks
                outputA = attackedOutputsA(i);
                attackList(outputA) = true;
            end
            obj.attackList = attackList;
        end
    end
end