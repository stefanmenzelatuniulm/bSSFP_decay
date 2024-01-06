classdef populationNode < emptyPopulationNode

    properties (Access=protected)

        parent populationNode;
    
    end

    methods

        %Constructor
        function populationNode = populationNode(parent)
            
            if nargin > 0
                populationNode.parent = parent;
            else
                populationNode.parent = emptyPopulationNode();
            end
    
        end

    end

end