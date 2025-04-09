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
        numSubObservers % number of sub observers per primary observer
        subsetIndices % which subobservers are a subobserver of a primary observer
        v0 % substation voltage function
        noise

        t % solution time
        x % solution states
    end

    methods
        function obj = PowerDistMO(numCustomers,attack,noise,v0,vref,LMIconsts,inputFileName)
            % PowerDistMO Construct an instance of this class
            %   Detailed explanation goes here

            fprintf(repmat('-',1,100));
            fprintf('\nCreating a system with %d customers\n',numCustomers);
            
            obj.numCustomers = numCustomers;
            obj.v0 = v0;
            obj.attack = attack;
            obj.noise = noise;

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

                assert(obj.attack.numAttacks < obj.numCustomers,...
                       "Number of attacks (%d) is larger than or equal to the number of customers (%d)",...
                       obj.attack.numAttacks,obj.numCustomers);
            end
            sysConsts.vref = vref;
            obj.sys = PowerSystem(numCustomers,sysConsts);

            numPrimaryObsvOutputs = obj.numCustomers - obj.attack.numAttacks;
            obj.primaryMO = MO(obj.sys,obj.attack,numPrimaryObsvOutputs,LMIconsts);
            numSecondaryObsvOutputs = 1;
            obj.secondaryMO = MO(obj.sys,obj.attack,numSecondaryObsvOutputs,LMIconsts);
            obj.numObservers = obj.primaryMO.numObservers + obj.secondaryMO.numObservers;

            [obj.numSubObservers,obj.subsetIndices] = obj.findIndices();

        end

        function [numOfSubObservers,subsetIndices] = findIndices(obj)
            % Find which secondary observers are a sub observer of all
            % primary observers.
            for j = 1:1:obj.primaryMO.numObservers
                CPrimIndices = obj.primaryMO.CSetIndices(j,:);
                % create new emtpy row to fill and append to the bottom of
                % PsubsetOfJIndices
                newRow = [];
                for p = 1:1:obj.secondaryMO.numObservers
                    CSecIndices = obj.secondaryMO.CSetIndices(p,:);
                    isPSubset = all(ismember(CSecIndices,CPrimIndices));
                    % If the indices of p are a subset of those of j: find the
                    if isPSubset
                        newRow(1,end+1) = p;
                    end
                end
                subsetIndices(j,:) = newRow;
            end
            numOfSubObservers = size(subsetIndices,2);
        end

        function [t,v,x,verr,xerr,bestObsv] = solve(obj,varargin) % tspan x0sys
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            tspan = [0 10];
            x0sys = zeros(1, obj.numCustomers);
            solveObservers = true;
            % extract values from varargin
            for i = 1:2:size(varargin,2)
                switch lower(varargin{i})
                    case 'tspan'
                        tspan = varargin{i+1};
                    case 'x0sys'
                        x0sys = varargin{i+1};
                    case 'solveobservers'
                        solveObservers = varargin{i+1};
                end
            end
            
            wb = waitbar(0,'Solver is currently at time: 0','Name','Solving the ODE');
            switch solveObservers
                case true
                    % construct x as [v; x; x_hat_primary; x_hat_secondary]
                    x0 = zeros(obj.sys.nx,1,1+obj.numObservers);
                    x0(:,:,1) = x0sys;
                    x0 = x0(:);

                    [obj.t,obj.x] = ode45(@(t,x) obj.odefunMO(wb,t,x,tspan(2)),tspan,x0);
                    % some reshape black magic to get it in the correct form (all
                    % observers behind each other)
                    obj.x = reshape(obj.x,[],obj.sys.nx,obj.numObservers+1); 
                    obj.x = permute(obj.x,[2 1 3]);
                    obj.t = obj.t';
                    t = obj.t;
                    xSys  = obj.x(:,:,1);

                    % select best estimate
                    delete(wb);
                    wb = waitbar(0,'Selection is currently at time: 0','Name','Selecting best estimates MO');
                    [bestObsv,xBestEst] = obj.selectBestEstimates(obj.t,obj.x,wb);
                    x = cat(3,xSys,xBestEst);
                    delete(wb);

                    numVoltageConversions = 2;

                case false
                    [obj.t,obj.x] = ode45(@(t,x) obj.odefunSys(wb,t,x,tspan(2)),tspan,x0sys);
                    % some reshape black magic to get it in the correct form (all
                    % observers behind each other)
                    obj.x = obj.x';
                    obj.t = obj.t';
                    x = obj.x;
                    t = obj.t;
                    numVoltageConversions = 1;
                    bestObsv = [];
                    delete(wb);
            end
            
            % calculate voltages

            v = zeros(obj.numCustomers,size(obj.t,2),numVoltageConversions);
            for k = 1:1:numVoltageConversions
                for ts = 1:1:size(obj.t,2)
                    v0sqrd = obj.v0(obj.t(ts))^2;
                    for i = 1:obj.numCustomers
                        % Stuff for v_i
                        obari = obj.obar(i);
                        v(i,ts,k) = obj.sys.C(i,:)*x(:,ts,1) + v0sqrd - obari;
                    end
                end
            end
            
            switch solveObservers
                case 1
                    verr = v(:,:,1) - v(:,:,2);
                    xerr = x(:,:,1) - x(:,:,2);
                case 0
                    verr = [];
                    xerr = [];
            end
        end
        
        function dx = odefunMO(obj,wb,t,x,tmax)
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
            xPrim = x(:,:,2:obj.primaryMO.numObservers+1);
            xSec  = x(:,:,obj.primaryMO.numObservers+2:end);
            
            % calculate u
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
                u(i) = obj.sys.Vref(t)^2 - obj.v0(t)^2 + oi;
            end

            % calculate system evolution
            mSys = obj.sys.C*xSys + u; % + disturbance
            ySys = mSys + obj.attack.value(t) + obj.noise.value(t);
            dxSys = obj.sys.A*xSys + obj.sys.B*obj.Q(mSys,[-14899.4 0 0 14899.4],obj.sys.charInverters);
            
            % calculate primary observers evolution
            dxPrim = zeros(size(xPrim));
            for o = 1:obj.primaryMO.numObservers
                xhat = xPrim(:,:,o);
                dyo = obj.primaryMO.CSet(:,:,o)*xhat - ySys(obj.primaryMO.CSetIndices(o,:));
                mhat = obj.sys.C*xhat + u + obj.primaryMO.K(:,:,o)*dyo;
                dxPrim(:,:,o) = obj.sys.A*xhat + obj.sys.B*obj.Q(mhat,[-14899.4 0 0 14899.4],obj.sys.charInverters) + obj.primaryMO.L(:,:,o)*dyo;
            end

            % calculate secondary observers evolution
            dxSec = zeros(size(xSec));
            for o = 1:obj.secondaryMO.numObservers
                xhat = xSec(:,:,o);
                dyo = obj.secondaryMO.CSet(:,:,o)*xhat - ySys(obj.secondaryMO.CSetIndices(o,:));
                mhat = obj.sys.C*xhat + u + obj.secondaryMO.K(:,:,o)*dyo;
                dxSec(:,:,o) = obj.sys.A*xhat + obj.sys.B*obj.Q(mhat,[-14899.4 0 0 14899.4],obj.sys.charInverters) + obj.secondaryMO.L(:,:,o)*dyo;
            end
            
            dx = cat(3,dxSys,dxPrim,dxSec);
            dx = dx(:);
        end

        function dx = odefunSys(obj,wb,t,x,tmax)
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

            % calculate u
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
                u(i) = obj.sys.Vref(t)^2 - obj.v0(t)^2 + oi;
            end

            % calculate system evolution
            m  = obj.sys.C*x + u; % + disturbance
            dx = obj.sys.A*x + obj.sys.B*obj.Q(m,[-14899.4 0 0 14899.4],obj.sys.charInverters);
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

        function [bestEstObsv,bestStateEstimate] = selectBestEstimates(obj,t,x,wb)
            xPrim = x(:,:,2:obj.primaryMO.numObservers+1);
            xSec  = x(:,:,obj.primaryMO.numObservers+2:end);
            % PiJ stores all the maximum difference between a primary and
            % all its secondary observers
            PiJ = zeros(obj.primaryMO.numObservers,1);
            
            tsteps = size(t,2);
            bestStateEstimate = zeros(obj.sys.nx,tsteps);
            bestEstObsv = zeros(1,tsteps);

            for t=1:tsteps
                % update waitbar
                try
                    waitbar(t/tsteps,wb,sprintf('Selector is currently at timestep: %d/%d',t,tsteps))
                catch ME
                    switch ME.identifier
                        case 'MATLAB:waitbar:InvalidSecondInput'
                            error('User terminated the solver.')
                        otherwise
                            rethrow(ME)
                    end
                end

                for primObsv = 1:obj.primaryMO.numObservers
                    xPrim_ti = xPrim(:,t,primObsv);
                    difflist = zeros(obj.numSubObservers,1);
                    for subObsv = 1:obj.numSubObservers
                        subObsvID = obj.subsetIndices(primObsv,subObsv);
                        xSub_ti = xSec(:,t,subObsvID);
                        difflist(subObsv) = norm(xPrim_ti-xSub_ti);
                    end
                    PiJ(primObsv) = max(difflist);
                end
                bestObsv = find(PiJ==min(PiJ));
                bestObsv = bestObsv(1);
                bestEstObsv(t) = bestObsv;
                bestStateEstimate(:,t) = xPrim(:,t,bestObsv);
            end
        end
    end
end