classdef longitudinalPopulation < population

    methods

        %Constructor
        function longitudinalPopulation = longitudinalPopulation(dephasingDegree,amplitude)
            
            if nargin > 0
                longitudinalPopulation = population(dephasingDegree,amplitude);
            else
                longitudinalPopulation = population();
            end

        end

    end

end