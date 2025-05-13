classdef Noise
    % class to define noise
    properties
        values
        times
    end

    methods
        function obj = Noise(inputs)
            % define number of timesamples for noise
            numOutputs = inputs.numCustomers;
            tspan = inputs.tspan;
            var = inputs.noiseVar;
            sampleFreq = inputs.sampleFreq;

            steps = (inputs.tspan(2) - tspan(1))*sampleFreq;
            stepsize = tspan(2)/steps;
            obj.times = tspan(1):stepsize:tspan(2);
            steps = size(obj.times,2);
            % setup noise vector
            obj.values = sqrt(var).*randn(numOutputs,steps);
        end

        function interpval = value(obj,t)
            % interpval returns the correct subset of outputs at the time t
            % the values will be interpolated from the correct timevalue if
            % t does not match with one that is in the obj.val
            interpval = interp1(obj.times',obj.values',t)';
        end

        function plot(obj)
            plot(obj.times,obj.values);
        end

    end

end