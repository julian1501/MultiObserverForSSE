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
        function obj = Attack(inputs)
            % Attack Construct an instance of this class
            %   Detailed explanation goes here
            numOutputs = inputs.numCustomers;
            obj.attackFunc = inputs.attack.func;
            obj.numAttacks = inputs.attack.num;
            if ~ (numOutputs > 2* obj.numAttacks)
               warning('The number of outputs is not larger then twice the number of attacked outputs %3.0f <= %3.0f',numOutputs,numAttacks); 
            end

            % pick random outputs to attack if none have been specified
            obj.attackedOutputs = inputs.attack.attackedOutputs;
            if isempty(obj.attackedOutputs)
                obj.attackedOutputs = sort(randperm(numOutputs,obj.numAttacks));
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