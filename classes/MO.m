classdef MO
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Sys % the system that this mo observes
        Attack % the attack that is used on the system
        numOutputs % number of system outputs (equal to ny if there are no sensor duplicates)
        numObservers % number of observers in this multi observer
        numOutputsObserver % number of outputs that each observer has
        Ci % Output matrix of each observer i
        CiIndices % Matrix containing the origin of each Ci
    end

    methods
        function obj = MO(Sys,Attack,numOutputs)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Sys = Sys;
            obj.Attack = Attack;
            obj.numOutputs = numOutputs;
            
            
        end
    end
end