classdef populationNode < emptyNode

    properties (Access = public)
    
        parent emptyNode;
        dephasingDegree double; %fraction of TR
        amplitude sym; %function of T1, T2, T2p
        amplitudeWithoutT2p sym; %T2p will be considered in populationTree based on the final dephasing degree
        amplitudeDirectlyAfterPulse sym;
        amplitudeDirectlyAfterPulseWithoutT2p sym; 
        amplitudeLabel sym;
        label string;
        level int64;
        xpos double;
        ypos double;
    
    end
    
    methods

        %Constructor
        function populationNode = populationNode(parent, label, xpos, ypos, dephasingDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p)

            if nargin > 1
    
                populationNode.parent = parent;
                populationNode.label = label;
                populationNode.dephasingDegree = dephasingDegree;
                populationNode.amplitude = amplitude;
                populationNode.amplitudeLabel = amplitudeLabel;

                if isa(parent, "populationNode")
                    populationNode.level = parent.level+1;
                else
                    populationNode.level = 0;
                end

                populationNode.xpos = xpos;
                populationNode.ypos = ypos;

                populationNode.amplitudeDirectlyAfterPulse = amplitudeDirectlyAfterPulse;
                populationNode.amplitudeWithoutT2p = amplitudeWithoutT2p;
                populationNode.amplitudeDirectlyAfterPulseWithoutT2p = amplitudeDirectlyAfterPulseWithoutT2p;

            else
                
                populationNode.parent = emptyNode();
                populationNode.label = "";
                populationNode.level = 0;
                populationNode.xpos = 0;
                populationNode.ypos = 0;
                populationNode.dephasingDegree = 0;
                populationNode.amplitude = sym(1);
                populationNode.amplitudeLabel = sym(1);
                populationNode.amplitudeDirectlyAfterPulse = sym(1);
                populationNode.amplitudeWithoutT2p = sym(1);
                populationNode.amplitudeDirectlyAfterPulseWithoutT2p = sym(1);

            end
    
        end

    end

end