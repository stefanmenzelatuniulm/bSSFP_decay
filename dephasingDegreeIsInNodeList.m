function isInNodeList = dephasingDegreeIsInNodeList(nodeList, dephasingDegree, level)

    isInNodeList = false;

    for k = 1:length(nodeList)
    
        if isa(nodeList(k), "populationNode")

            populationNode = nodeList(k);

            if populationNode.level == level

                if ismembertol(populationNode.dephasingDegree, dephasingDegree, 0.001)
    
                    isInNodeList = true;
                    break;
    
                end

            end

        end
    
    end

end