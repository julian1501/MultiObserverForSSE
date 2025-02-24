classdef PowerSystem
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        A % State matrix
        B % Input matrix
        C % Output matrix
        nx
        ny
        nu
    end

    methods
        function obj = PowerSystem(numCustomers, inverters, activeImpedances, reactiveImpedances)
            % PowerSystem Construct an instance of this class
            obj.A = diag(-1.*inverters);
            obj.B = diag(inverters);

            activeImpMainLine   =   activeImpedances(:,1);
            activeImpBranches   =   activeImpedances(:,2);
            
            % Generate C matrix template
            obj.C = zeros(numCustomers,numCustomers);
            for i = 1:1:numCustomers
                obj.C(i:end,i:end) = obj.C(i:end,i:end) + repmat(activeImpMainLine(i),numCustomers+1-i,numCustomers+1-i);
            end
            obj.C = obj.C - diag(2.*activeImpBranches);

            obj.nx = size(A,1);
            obj.ny = size(C,1);
            obj.nu = size(B,2);
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