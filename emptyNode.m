classdef emptyNode

    %Important for recursion exit

    methods (Static)

        function [plottedTransverseNodes, plottedLongitudinalNodes] = plotNode(~, ~, ~, ~, ~, ~, ~, ~, ~)

            plottedTransverseNodes = [];
            plottedLongitudinalNodes = [];

        end

        function plotPathway(~, ~)
        end

    end

    methods

        function emptyNodeObject = prune(emptyNodeObject, ~)            
        end

        function emptyNodeObject = updateAmplitudeLabel(emptyNodeObject, ~, ~, ~, ~, ~, ~, ~)
        end

    end

end