classdef populationNode < emptyNode

    properties (Access=public)
    
        parent emptyNode;
        dephasingDegree double; %fraction of TR
        amplitude sym; %function of T1, T2, T2p
        label string;
        level int64;
        xpos double;
        ypos double;
        Meq double;
    
    end
    
    methods

        %Constructor
        function populationNode = populationNode(parent,label,xpos,ypos,Meq,dephasingDegree,amplitude)

            if nargin > 1
    
                populationNode.parent = parent;
                populationNode.label = label;
                populationNode.dephasingDegree = dephasingDegree;
                populationNode.amplitude = amplitude;

                if isa(parent,"populationNode")
                    populationNode.level = parent.level+1;
                else
                    populationNode.level = 0;
                end

                populationNode.xpos = xpos;
                populationNode.ypos = ypos;
                populationNode.Meq = Meq;

            else
                
                populationNode.parent = emptyNode();
                populationNode.label = "";
                populationNode.level = 0;
                populationNode.xpos = 0;
                populationNode.ypos = 0;
                populationNode.Meq = 1;
                populationNode.dephasingDegree = 0;
                populationNode.amplitude = sym(1);

            end
    
        end

    end

end