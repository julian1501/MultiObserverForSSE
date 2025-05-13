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
        vref % reference voltage
        noise
        predictors % boolean indicating whether predictors are active or not
        interSampleTimes % sample intervals between sensor updates
    end

    methods
        function obj = PowerDistMO(inputs,attack,noise)
            % PowerDistMO Construct an instance of this class
            %   Detailed explanation goes here

            fprintf(repmat('-',1,100));
            fprintf('\nCreating a system with %d customers\n',obj.numCustomers);
            
            obj.numCustomers = inputs.numCustomers;
            obj.v0 = inputs.sys.v0;
            obj.attack = attack;
            obj.noise = noise;
            obj.predictors = inputs.predictors;

            if ~isempty(inputs.inputFileName)
                ... % read file extract data
                % Implement when data format is known
                error('File input not yet supported')
            else
                % in case of no file provided use 'standard system'
                sysConsts = struct();
                sysConsts.inverters = ones(obj.numCustomers,1);
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

            obj.vref = inputs.sys.vref;
            sysConsts.vref = inputs.sys.vref;
            obj.sys = PowerSystem(obj.numCustomers,sysConsts);
            obsvSize = inputs.obsv.size;

            if obsvSize(1) == 0
                primObsvOutputs = obj.numCustomers - obj.attack.numAttacks;
            else
                primObsvOutputs = obsvSize(1);
            end
            if obsvSize(2) == 0
                if attack.numAttacks*2 < obj.numCustomers
                    secObsvOutputs = obj.numCustomers - 2*obj.attack.numAttacks;
                else
                    secObsvOutputs  = 1;
                end
            else
                secObsvOutputs = obsvSize(2);
            end
            
            obj.primaryMO = MO(obj.sys,obj.attack,primObsvOutputs,inputs.LMIconsts);
            obj.secondaryMO = MO(obj.sys,obj.attack,secObsvOutputs,inputs.LMIconsts);
            obj.numObservers = obj.primaryMO.numObservers + obj.secondaryMO.numObservers;

            [obj.numSubObservers,obj.subsetIndices] = obj.findIndices();

            obj.interSampleTimes = obj.initPredictors(5);

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

        function outputs = solve(obj,varargin) % tspan x0sys
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
            
            if obj.predictors.enabled == 1
                E = odeEvent('EventFcn',@obj.predictorUpdateCheck, ...
                             'Direction',"descending", ...
                             'Response',"callback", ...
                             'CallbackFcn',@updateWrapper);
            else
                E = [];
            end

            % construct x as [v; x; x_hat_primary; x_hat_secondary]
            x0 = zeros(obj.sys.nx,1,3+obj.numObservers);
            x0(:,:,1) = obj.interSampleTimes(:,2);
            x0(:,:,2) = x0sys;

            % setup ode
            odeFun = @(t,x,p) obj.odefunMO(t,x,p);
            ODE = ode(ODEFcn=odeFun,InitialValue=x0(:),EventDefinition=E,Solver="ode45");
            ODE.Parameters = {wb, tspan(2), obj};
            sol = solve(ODE,tspan(1),tspan(2));
            outputs.t = sol.Time;
            xs = sol.Solution;

            % get the observers behind each other
            x = zeros(obj.numCustomers,size(outputs.t,2),obj.numObservers + 3);
            for i = 1:1:obj.numObservers + 3
                rows = (i-1)*obj.numCustomers + (1:obj.numCustomers);
                x(:,:,i) = xs(rows,:);
            end
            % x = permute(x,[2 1 3]);
            % t = t';             
            
            outputs.timers = x(:,:,1);
            % recover y from x
            outputs.yhat = x(:,:,3);
            outputs.y = zeros(size(x(:,:,2)));
            for ts = 1:1:size(outputs.t,2)
                outputs.y(:,ts) = obj.sys.C*x(:,ts,2) + obj.u(outputs.t(ts)) + ...
                    obj.attack.value(outputs.t(ts)) + obj.noise.value(outputs.t(ts));
            end
           
            % select best estimate
            delete(wb);
            wb = waitbar(0,'Selection is currently at time: 0','Name','Selecting best estimates MO');
            [outputs.bestObsv,outputs.xBestEst] = obj.selectBestEstimates(outputs.t,x,wb);
            outputs.xout = cat(3,x(:,:,2),outputs.xBestEst);
            delete(wb);

            numVoltageConversions = 2;
            
            % calculate voltages

            outputs.vout = zeros(obj.numCustomers,size(outputs.t,2),numVoltageConversions);
            for k = 1:1:numVoltageConversions
                for ts = 1:1:size(outputs.t,2)
                    v0sqrd = obj.v0(outputs.t(ts))^2;
                    for i = 1:obj.numCustomers
                        % Stuff for v_i
                        obari = obj.obar(i);
                        outputs.vout(i,ts,k) = obj.sys.C(i,:)*outputs.xout(:,ts,k) + v0sqrd - obari;
                    end
                end
            end
            
            switch solveObservers
                case 1
                    outputs.verr = outputs.vout(:,:,1) - outputs.vout(:,:,2);
                    outputs.xerr = outputs.xout(:,:,1) - outputs.xout(:,:,2);
                case 0
                    outputs.verr = [];
                    outputs.xerr = [];
            end
        end
        
        function dx = odefunMO(obj,t,x,p)
            wb = p{1,1};
            tmax = p{1,2};
            obj = p{1,3};

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
            timer = x(:,:,1);
            xSys  = x(:,:,2);
            yhat  = x(:,:,3);
            xPrim = x(:,:,4:obj.primaryMO.numObservers+3);
            xSec  = x(:,:,obj.primaryMO.numObservers+4:end);
            
            u = obj.u(t);

            % set timer derivatives
            dtimer = -1.*ones(size(timer));

            % calculate system evolution
            mSys = obj.sys.C*xSys + u;
            dxSys = obj.sys.A*xSys + obj.sys.B*obj.Q(mSys,[-14899.4 0 0 14899.4],obj.sys.charInverters);
            ySys = mSys + obj.attack.value(t) + obj.noise.value(t);
            
            % update yhat
            [~,xhat] = obj.selectBestEstimates(1,x,0);
            mSysPred = obj.sys.C*xhat + u;
            Bphi = obj.sys.B*obj.Q(mSysPred,[-14899.4 0 0 14899.4],obj.sys.charInverters);
            dyhat = obj.sys.C*(obj.sys.A*xhat + Bphi) - obj.predictors.gain*eye(obj.numCustomers)*(yhat - obj.sys.C*xhat);

            
            % calculate primary observers evolution
            dxPrim = zeros(size(xPrim));
            for o = 1:obj.primaryMO.numObservers
                xhat = xPrim(:,:,o);
                dyo  = obj.primaryMO.CSet(:,:,o)*xhat + u(obj.primaryMO.CSetIndices(o,:)) - yhat(obj.primaryMO.CSetIndices(o,:));
                mhat = obj.sys.C*xhat + u + obj.primaryMO.K(:,:,o)*dyo;
                dxPrim(:,:,o) = obj.sys.A*xhat + obj.sys.B*obj.Q(mhat,[-14899.4 0 0 14899.4],obj.sys.charInverters) + obj.primaryMO.L(:,:,o)*dyo;
            end

            % calculate secondary observers evolution
            dxSec = zeros(size(xSec));
            for o = 1:obj.secondaryMO.numObservers
                xhat = xSec(:,:,o);
                dyo  = obj.secondaryMO.CSet(:,:,o)*xhat + u(obj.secondaryMO.CSetIndices(o,:)) - yhat(obj.secondaryMO.CSetIndices(o,:));
                mhat = obj.sys.C*xhat + u + obj.secondaryMO.K(:,:,o)*dyo;
                dxSec(:,:,o) = obj.sys.A*xhat + obj.sys.B*obj.Q(mhat,[-14899.4 0 0 14899.4],obj.sys.charInverters) + obj.secondaryMO.L(:,:,o)*dyo;
            end
            
            dx = cat(3,dtimer,dxSys,dyhat,dxPrim,dxSec);
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
            u = obj.u(t);

            % calculate system evolution
            m  = obj.sys.C*x + u; % + disturbance
            dx = obj.sys.A*x + obj.sys.B*obj.Q(m,[-14899.4 0 0 14899.4],obj.sys.charInverters);
        end

        function minTimer = predictorUpdateCheck(obj,t,x)
            x = reshape(x,obj.sys.nx,1,[]);
            timers = x(:,:,1);            
            minTimer = min(timers);
            sensorToUpdate = find(timers==minTimer);
            
            if minTimer < 0
                fprintf("Sensor %d requires update at %2.2f [s]\n",sensorToUpdate,t)
            end
        end

        function [stop, xupdated] = predictorUpdate(obj,t,x,i,parameters)
            stop = false; % do not stop the simulation

            x = reshape(x,obj.sys.nx,1,[]);
            timers = x(:,:,1);
            xSys   = x(:,:,2);
            yhat   = x(:,:,3);
            y = obj.sys.C*xSys + obj.u(t) + obj.attack.value(t) + obj.noise.value(t);

            minTimer = min(timers);
            sensorToUpdate = find(timers==minTimer);

            yhat(sensorToUpdate) = y(sensorToUpdate);
            timers(sensorToUpdate) = obj.interSampleTimes(sensorToUpdate,1);
            
            xupdated = x;
            xupdated(:,:,1) = timers;
            xupdated(:,:,3) = yhat;
            xupdated = xupdated(:);

            fprintf("Sensor %d updated at %2.2f\n",sensorToUpdate,t)
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

        function u = u(obj,t)
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
        end

        function [bestEstObsv,bestStateEstimate] = selectBestEstimates(obj,t,x,wb)
            xPrim = x(:,:,3:obj.primaryMO.numObservers+2);
            xSec  = x(:,:,obj.primaryMO.numObservers+3:end);
            % PiJ stores all the maximum difference between a primary and
            % all its secondary observers
            tsteps = max(size(t));
            bestStateEstimate = zeros(obj.sys.nx,tsteps);
            bestEstObsv = zeros(1,tsteps);

            for ts=1:tsteps
                % update waitbar
                if wb ~= 0
                    try
                        waitbar(ts/tsteps,wb,sprintf('Selector is currently at timestep: %d/%d',ts,tsteps))
                    catch ME
                        switch ME.identifier
                            case 'MATLAB:waitbar:InvalidSecondInput'
                                error('User terminated the solver.')
                            otherwise
                                rethrow(ME)
                        end
                    end
                end
                
                PiJ = zeros(obj.primaryMO.numObservers,1);
                for primObsv = 1:obj.primaryMO.numObservers
                    xPrim_ti = xPrim(:,ts,primObsv);
                    difflist = zeros(obj.numSubObservers,1);
                    for subObsv = 1:obj.numSubObservers
                        subObsvID = obj.subsetIndices(primObsv,subObsv);
                        xSub_ti = xSec(:,ts,subObsvID);
                        difflist(subObsv) = norm(xPrim_ti-xSub_ti);
                    end
                    PiJ(primObsv) = max(difflist);
                end
                bestObsv = find(PiJ==min(PiJ));
                bestObsv = bestObsv(1);
                bestEstObsv(ts) = bestObsv;
                bestStateEstimate(:,ts) = xPrim(:,ts,bestObsv);
            end
        end

        function interSampleTimes = initPredictors(obj,maxTinterval)
            
            % column 1: inter sample time
            % column 2: start offset (< inter sample time)
            interSampleTimes = zeros(obj.numCustomers,2); 
            for c = 1:1:obj.numCustomers
                interSampleTimes(c,1) = rand(1)*maxTinterval;
                interSampleTimes(c,2) = rand(1)*interSampleTimes(c,1);
            end
            % normalize times such that the first update is at t=0
            minStartTime = min(interSampleTimes(:,2));
            % add offset to prevent ode event from immediately triggering
            interSampleTimes(:,2) = interSampleTimes(:,2) - minStartTime + 0.01; 
        end
    end
end

function [stop, xupdated] = updateWrapper(t,x,i,p)
    obj = p{1,3};
    [stop, xupdated] = obj.predictorUpdate(t,x,[],p);
end
