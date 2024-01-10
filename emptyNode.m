classdef emptyNode

    %Important for recursion exit

    methods (Static)

        function [plottedTransverseDephasingDegrees, plottedLongitudinalDephasingDegrees] = plotNode(~, ~, ~, ~, ~, ~, ~, ~, ~, ~)

            plottedTransverseDephasingDegrees = [];
            plottedLongitudinalDephasingDegrees = [];

        end

        function plotPathway(~, ~)
        end

    end

    methods

        function emptyNodeObject = prune(emptyNodeObject, ~)            
        end

        function emptyNodeObject = updateAmplitudeLabel(emptyNodeObject, ~, ~, ~, ~)
        end

    end

end