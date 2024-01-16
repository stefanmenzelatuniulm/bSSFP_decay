function isInNodeList = coherenceDegreeIsInNodeList(nodeList, coherenceDegree, level)

    isInNodeList = false;

    for k = 1:length(nodeList)
    
        if isa(nodeList(k), "populationNode")

            populationNode = nodeList(k);

            if populationNode.level == level

                if ismembertol(populationNode.coherenceDegree, coherenceDegree, 0.001)
    
                    isInNodeList = true;
                    break;
    
                end

            end

        end
    
    end

end