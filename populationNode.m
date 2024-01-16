classdef populationNode < emptyNode

    properties (Access = public)
    
        parent emptyNode;
        coherenceDegree double; 
        coherenceDegreeDirectlyAfterPulse double;
        amplitude sym; %function of T1, T2, T2p
        amplitudeWithoutT2p sym; %T2p will be considered in populationTree based on the final dephasing degree
        amplitudeDirectlyAfterPulse sym;
        amplitudeDirectlyAfterPulseWithoutT2p sym; 
        amplitudeLabel sym;
        label string;
        level double;
        totalTime double;
    
    end
    
    methods

        %Constructor
        function populationNode = populationNode(parent, label, totalTime, coherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulse)

            if nargin > 1
    
                populationNode.parent = parent;
                populationNode.label = label;
                populationNode.coherenceDegree = coherenceDegree;
                populationNode.amplitude = amplitude;
                populationNode.amplitudeLabel = amplitudeLabel;

                if isa(parent, "populationNode")
                    populationNode.level = parent.level+1;
                else
                    populationNode.level = 0;
                end

                populationNode.totalTime = totalTime;

                populationNode.amplitudeDirectlyAfterPulse = amplitudeDirectlyAfterPulse;
                populationNode.amplitudeWithoutT2p = amplitudeWithoutT2p;
                populationNode.amplitudeDirectlyAfterPulseWithoutT2p = amplitudeDirectlyAfterPulseWithoutT2p;
                populationNode.coherenceDegreeDirectlyAfterPulse = coherenceDegreeDirectlyAfterPulse;

            else
                
                populationNode.parent = emptyNode();
                populationNode.label = "";
                populationNode.level = 0;
                populationNode.totalTime = 0;
                populationNode.coherenceDegree = 0;
                populationNode.amplitude = sym(1);
                populationNode.amplitudeLabel = sym(1);
                populationNode.amplitudeDirectlyAfterPulse = sym(1);
                populationNode.amplitudeWithoutT2p = sym(1);
                populationNode.amplitudeDirectlyAfterPulseWithoutT2p = sym(1);
                populationNode.coherenceDegreeDirectlyAfterPulse = 0;

            end
    
        end

    end

end