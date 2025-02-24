classdef MO
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        sys % the system that this mo observes
        attack % the attack that is used on the system
        numOutputs % the number of sensors sensing the system
        numObservers % number of observers in this multi observer
        numOutputsObserver % number of outputs that each observer has
        CSet % Output matrix of each observer i
        CSetIndices % Matrix containing the origin of each Ci
    end

    methods
        function obj = MO(sys,attack)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
            obj.attack = attack;
            obj.numOutputs = obj.sys.nx;

            if 2*obj.attack.numAttacks > obj.numOutputs
                warning(["More than half of the sensors are attacked",...
                        "2*numAttacks (%d) > numOutputs (%d)"],...
                        obj.attack.numAttacks,obj.numOutputs);
            end
            
            obj.numOutputsObserver  = obj.numOutputs - obj.attack.numAttacks;
            obj.numObservers        = nchoosek(obj.numOutputs,obj.numOutputsObserver);
            
            [obj.CSet,obj.CSetIndices] = CSetSetup(obj);

        end

        function [CSet,CSetIndices] = CSetSetup(obj)
            
            CSetIndices = nchoosek(1:obj.numOutputs,obj.numOutputsObserver);
            CSet = zeros(obj.numOutputsObserver,...
                         obj.sys.nx,...
                         obj.numObservers);
            for i = 1:obj.numObservers
                CSet(:,:,i) = obj.sys.C(CSetIndices(i,:),:);
            end
        end
    end
end