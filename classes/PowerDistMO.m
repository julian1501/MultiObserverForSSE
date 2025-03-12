classdef PowerDistMO
    % PowerDistMO Summary of this class goes here
    %   Detailed explanation goes here

    properties
        numCustomers % Number of customers in the system
        sys % Power Distribution system object
        attack % Attack object
        primaryMO % Large multi observer
        secondaryMO % Small multi observer
        numObservers % Total number of observers
        powerGenCon % power consumption/generation data
        sysConsts
        t % solution time
        x % solution states
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
                sysConsts = struct();
                sysConsts.inverters = ones(numCustomers,1);
                sysConsts.actImp = [0.00343, 0.04711; 
                                    0.00172, 0.02356; 
                                    0.00343, 0.04711; 
                                    0.00515, 0.07067; 
                                    0.00172, 0.02356];
                
                sysConsts.reaImp   = [0.00147, 0.02157;
                                      0.00662, 0.09707;
                                      0.00147, 0.02157;
                                      0.00147, 0.02157;
                                      0.00147, 0.02157];
                
                obj.powerGenCon.actGen = [3500 5500 4000 4500 3000];
                obj.powerGenCon.actCon = [2295 5440 5440 2295 2720];
                obj.powerGenCon.reaCon = [ 300  960  480  600  400];
                
                sysConsts.charInverters = [2321.6 3464.1 2467.8 2800.0 1900.0];

                if obj.numCustomers ~= 5
                    reps = ceil(obj.numCustomers/5);
                    sysConsts.inverters = repmat(sysConsts.inverters,reps,1);
                    sysConsts.inverters = sysConsts.inverters(1:obj.numCustomers);
                    sysConsts.actImp    = repmat(sysConsts.actImp,reps,1);
                    sysConsts.actImp    = sysConsts.actImp(1:obj.numCustomers,:);
                    sysConsts.reaImp    = repmat(sysConsts.reaImp,reps,1);
                    sysConsts.reaImp    = sysConsts.reaImp(1:obj.numCustomers,:);
                    obj.powerGenCon.actGen = repmat(obj.powerGenCon.actGen,1,reps);
                    obj.powerGenCon.actGen = obj.powerGenCon.actGen(:,1:obj.numCustomers);
                    obj.powerGenCon.actCon = repmat(obj.powerGenCon.actCon,1,reps);
                    obj.powerGenCon.actCon = obj.powerGenCon.actCon(:,1:obj.numCustomers);
                    obj.powerGenCon.reaCon = repmat(obj.powerGenCon.reaCon,1,reps);
                    obj.powerGenCon.reaCon = obj.powerGenCon.reaCon(:,1:obj.numCustomers);
                    sysConsts.charInverters = repmat(sysConsts.charInverters,1,reps);
                    sysConsts.charInverters = sysConsts.charInverters(:,1:obj.numCustomers);
                end

                assert(numAttacks < obj.numCustomers,...
                       "Number of attacks (%d) is larger than or equal to the number of customers (%d)",...
                       numAttacks,obj.numCustomers);
            end
            
            obj.sys = PowerSystem(numCustomers,sysConsts);
            obj.attack = Attack(obj.numCustomers,numAttacks,attackedOutputs);

            numPrimaryObsvOutputs = obj.attack.numAttacks;
            obj.primaryMO = MO(obj.sys,obj.attack,numPrimaryObsvOutputs);
            numSecondaryObsvOutputs = 1;
            obj.secondaryMO = MO(obj.sys,obj.attack,numSecondaryObsvOutputs);
            obj.numObservers = obj.primaryMO.numObservers + obj.secondaryMO.numObservers;

            [numOfSubObservers,subsetIndices] = obj.findIndices();

        end

        function [numOfSubObservers,subsetIndices] = findIndices(obj)
            % Find which secondary observers are a sub observer of all
            % primary observers.
            for j = 1:1:obj.primaryMO.numObservers
                CjIndices = obj.primaryMO.CiIndices(j,:);
                % create new emtpy row to fill and append to the bottom of
                % PsubsetOfJIndices
                newRow = [];
                for p = 1:1:obj.secondaryMO.numObservers
                    CpIndices = Pmo.CiIndices(p,:);
                    isPSubset = isSubsetOf(CjIndices,CpIndices);
                    % If the indices of p are a subset of those of j: find the
                    if isPSubset
                        newRow(1,end+1) = p;
                    end
                end
                subsetIndices(j,:) = newRow;
            end
            numOfSubObservers = size(PsubsetOfJIndices,2);
        end

        function [t,v,x] = solve(obj,tspan,x0sys)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            % construct x as [v; x; x_hat_primary; x_hat_secondary]
            x0 = zeros(obj.sys.nx,1,1+obj.numObservers);
            x0(:,:,1) = x0sys;
            x0 = x0(:,:,1);
            x0 = x0(:);

            wb = waitbar(0,'Solver is currently at time: 0','Name','Solving the ODE');
            [t,x] = ode45(@(t,x) obj.odefun(wb,t,x,tspan(2)),tspan,x0);
            x = x'; t = t';
            xSys = x(1:obj.numCustomers,:);
            % extract v_is
            v = zeros(obj.numCustomers,size(t,2));
            for ts = 1:1:size(t,2)
                for i = 1:obj.numCustomers
                    % Stuff for v_i
                    obari = obj.obar(i);
                    v(i,ts) = sqrt(obj.sys.C(i,:)*xSys(:,ts) + obj.v0(t(ts))^2 - obari);
                end
            end
            delete(wb);
        end
        
        function dx = odefun(obj,wb,t,x,tmax)
            % Update waitbar
            try
                waitbar(t/tmax,wb,sprintf('Solver is currently at time: %.4f',t))
            catch ME
                switch ME.identifier
                    case 'MATLAB:waitbar:InvalidSecondInput'
                        error('User terminated the solver.')
                    otherwise
                        rethrow(ME)
                end
        
            end
            % Extract different x sets
            x = reshape(x,obj.sys.nx,1,[]);
            xSys  = x(:,:,1);

%             xPrim = x(:,:,2:obj.secondaryMO.numObservers+1);
%             xSec  = x(:,:,obj.secondaryMO.numObservers+2:end);
            
            % calculate system evolution
            u = zeros(obj.numCustomers,1);
            for i = 1:obj.numCustomers
                oi = 0;

                % Stuff for u_i
                for j = 1:i
                    obarj = obj.obar(j);
                    if ~(j < i)
                        oi = oi + obarj;
                        continue;
                    end
                    betaj = obj.beta(j);
                
                    oi = oi + obarj + 2*betaj;
                end
                u(i) = obj.sys.Vref^2 - obj.v0(t)^2 + oi;
            end
            mSys = obj.sys.C*xSys + u; % + d
            dxSys = obj.sys.A*xSys + obj.sys.B*obj.Q(mSys,[-14899.4 0 0 14899.4],obj.sys.charInverters);
            dx = dxSys;
        end

        function Q_m = Q(obj,m,wLims,refCon)
            % Q function: Saturated function with a deadzone
            % m: Input vector (numCustomers x 1)
            % mlims: Saturation and deadzone limits (4 x 1)
            % refCon: reference constants
            % Q_m: Output vector (numCustomers x 1)
            Q_m = zeros(obj.numCustomers,1);
            for i = 1:1:obj.numCustomers
                if m(i) <= wLims(1)
                    Q_m(i) = -refCon(i);
                elseif (wLims(1) < m(i)) && (m(i) <= wLims(2))
                    Q_m(i) = -(1-(m(i)-wLims(1))/(wLims(2)-wLims(1)))*refCon(i);
                elseif (wLims(2) < m(i)) && (m(i) <= wLims(3))
                    Q_m(i) = 0;
                elseif (wLims(3) < m(i)) && (m(i) <= wLims(4))
                    Q_m(i) = ((m(i)-wLims(3))/(wLims(4)-wLims(3)))*refCon(i);
                elseif wLims(4) <= m(i)
                    Q_m(i) = refCon(i);
                else
                    error('Nonlinearity failed')
                end
            end
        end

        function v0t = v0(~,t)
            v0t = 230 + 5*sin(t);
        end

        function obari = obar(obj,i)
            % Implement obar function as in the paper
            reactivePower = 0;
            activePower   = 0;
            for k = i:obj.numCustomers
                reactivePower = reactivePower + obj.powerGenCon.reaCon(k);
                activePower   = activePower   + obj.powerGenCon.actGen(k) - obj.powerGenCon.actCon(k);
            end
            reactive = obj.sys.reaImp(i,1)*reactivePower;
            active   = obj.sys.actImp(i,1)*activePower;
            obari = 2*reactive - 2*active - 2*obj.beta(i);

        end

        function betai = beta(obj,i)

            if i == 0
                betai = 0;
            elseif i > 0
                betai = obj.sys.actImp(i,2)*(obj.powerGenCon.actGen(i) - obj.powerGenCon.actCon(i)) + ...
                        obj.sys.reaImp(i,2)*obj.powerGenCon.reaCon(i);
            end
        end
    end
end