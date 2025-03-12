classdef Attack
    % Attack Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numAttacks
        attackList
        attackFunc
        attackedOutputs
    end

    methods
        function obj = Attack(numOutputs,numAttacks,attackedOutputs,attackFunc)
            % Attack Construct an instance of this class
            %   Detailed explanation goes here
            obj.attackFunc = attackFunc;
            obj.numAttacks = numAttacks;
            if ~ (numOutputs > 2* numAttacks)
               warning('The number of outputs is not larger then twice the number of attacked outputs %3.0f <= %3.0f',numOutputs,numAttacks); 
            end

            % pick random outputs to attack if none have been specified
            if isempty(attackedOutputs)
                obj.attackedOutputs = sort(randperm(numOutputs,numAttacks));
            end
            fprintf("The attacked outputs are: \n");
            disp(obj.attackedOutputs);
            
            % loop over outputs and set booleans indicating if it is
            % attacked
            obj.attackList = false(numOutputs,1);
            obj.attackList(obj.attackedOutputs) = true;
            
        end

        function attackVal = value(obj,t)
            % Calculate the value of the attack signal for each output that
            % is attacked.
            attackVal = zeros(size(obj.attackList));
            attackVal(obj.attackList) = obj.attackFunc(t);
        end
    end
end