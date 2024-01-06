classdef population

    properties
        
        dephasingDegree double; %fraction of TR
        amplitude sym; %function of T1, T2, T2p

    end

    methods

        %Constructor
        function population = population(dephasingDegree,amplitude)
            
            if nargin > 0
                population.dephasingDegree = dephasingDegree;
                population.amplitude = amplitude;
            else
                population.dephasingDegree = 0;
                population.amplitude = 0;
            end

        end

    end

end