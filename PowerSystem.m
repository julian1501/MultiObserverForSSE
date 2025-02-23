classdef PowerSystem
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A % State matrix
        B % Input matrix
        C % Output matrix
    end

    methods
        function obj = PowerSystem(numCustomers, inverters, activeImpedances, reactiveImpedances)
            % PowerSystem Construct an instance of this class
            %   Detailed explanation goes here
            obj.A = diag(-1.*inverters);
            obj.B = diag(inverters);
            
            % Generate C matrix template
            Ctemplate = zeros(numCustomers,numCustomers);
            for i = 1:1:numCustomers
                Ctemplate(i:end,i:end) = Ctemplate(i:end,i:end) + repmat(X(i),numCustomers-1,numCustomers-1);
            end

            % Generate C matrix by slicing the templates
            C = zeros(numCustomers,numCustomers,numCustomers);
            for i = 1:1:numCustomers
                C(1:i,1:i,i) = Ctemplate(1:i,1:i);
            end
        end

        function stableBool = isStable(obj)
            % isStable
            % Check if system is stable
            eigenvalues = eig(obj.A);
            if any(eigenvalues > 0)
                stableBool = false;
            else
                stableBool = true;
            end
        end
    end
end