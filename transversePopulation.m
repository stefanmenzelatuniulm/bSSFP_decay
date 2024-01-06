classdef transversePopulation < population

    methods

        %Constructor
        function transversePopulation = transversePopulation(dephasingDegree,amplitude)
            
            if nargin > 0
                transversePopulation = population(dephasingDegree,amplitude);
            else
                transversePopulation = population();
            end

        end

    end

end