classdef Noise
    % class to define noise
    properties
        values
        times
    end

    methods
        function obj = Noise(numOutputs,tspan,var,sampleFreq)
            % define number of timesamples for noise
            steps = (tspan(2) - tspan(1))*sampleFreq;
            stepsize = tspan(2)/steps;
            obj.times = tspan(1):stepsize:tspan(2);
            steps = size(obj.times,2);
            % setup noise vector
            obj.values = sqrt(var).*randn(numOutputs,steps);
        end

        function interpval = interpNoise(obj,outputs,t)
            % interpval returns the correct subset of outputs at the time t
            % the values will be interpolated from the correct timevalue if
            % t does not match with one that is in the obj.val
            %
            % if outputs is set to the string 'all' all outputs are used
            if strcmp(outputs,"all")
                outputs = 1:1:size(obj.values,1);
            end
            interpFullData = interp1(obj.times',obj.values',t)';
            interpval = interpFullData(outputs,1);
        end

    end

end