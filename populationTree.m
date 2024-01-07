classdef populationTree

    properties (Access = public)

        root longitudinalPopulationNode;
        height int64;
        Meq double;
        yScale double;

    end
        
    methods
        
        %Constructor
        function populationTree = populationTree(root, Meq, yScale)

            if nargin > 1
            
                populationTree.root = root;
                populationTree.Meq = Meq;
                populationTree.yScale = yScale;

            else

                populationTree.root = emptyNode();
                populationTree.Meq = 1;
                populationTree.yScale = 1;

            end

            populationTree.height = 0;
    
        end

        %Applies pulse to population tree
        function [transverseBottomNodes, longitudinalBottomNodes, populationTree] = applyPulse(populationTree, a, TR, f)
            
            [transverseBottomNodes, longitudinalBottomNodes, populationTree.root] = populationTree.root.applyPulse(a, TR, f, populationTree.height, populationTree.Meq, populationTree.yScale); %bottomNodes is list with entries label+"#"+string(amplitude)+"#"+string(dephasingDegree)
            populationTree.height = populationTree.height+1;

        end

        %Sums branches with the same dephasing degree by assigning the
        %summed amplitudes to one of the branches, pruning the rest of
        %the branches and updating labels
        function populationTree = pruneMerge(populationTree, transverseBottomNodes, longitudinalBottomNodes)

            transverseBottomNodes = transpose(split(transpose(transverseBottomNodes), "#"));
            longitudinalBottomNodes = transpose(split(transpose(longitudinalBottomNodes), "#"));
            
            transverseLabels = transverseBottomNodes(:, 1);
            transverseAmplitudes = str2sym(transverseBottomNodes(:, 2));
            transverseDephasingDegrees = double(transverseBottomNodes(:, 3));

            longitudinalLabels = longitudinalBottomNodes(:, 1);
            longitudinalAmplitudes = str2sym(longitudinalBottomNodes(:, 2));
            longitudinalDephasingDegrees = double(longitudinalBottomNodes(:, 3));
            
            [longitudinalUpdateIndices, longitudinalPruneIndices, longitudinalSummedAmplitudes, longitudinalUpdateFullLabels] = populationTreeMergePruneHelper(longitudinalDephasingDegrees, longitudinalAmplitudes, longitudinalLabels);
            longitudinalPruneLabels = longitudinalLabels(longitudinalPruneIndices);
            longitudinalUpdateLabels = longitudinalLabels(longitudinalUpdateIndices);

            [transverseUpdateIndices, transversePruneIndices, transverseSummedAmplitudes, transverseUpdateFullLabels] = populationTreeMergePruneHelper(transverseDephasingDegrees, transverseAmplitudes, transverseLabels);
            transversePruneLabels = transverseLabels(transversePruneIndices);
            transverseUpdateLabels = transverseLabels(transverseUpdateIndices);

            for k = 1:length(longitudinalPruneLabels)
            
                populationTree = populationTree.prune(longitudinalPruneLabels(k));
            
            end

            for k = 1:length(longitudinalUpdateLabels)
                populationTree = populationTree.updateAmplitudeLabel(longitudinalUpdateLabels(k), longitudinalSummedAmplitudes(k), longitudinalUpdateFullLabels(k));
            
            end

            for k = 1:length(transversePruneLabels)
            
                populationTree = populationTree.prune(transversePruneLabels(k));
            
            end

            for k = 1:length(transverseUpdateLabels)
            
                populationTree = populationTree.updateAmplitudeLabel(transverseUpdateLabels(k), transverseSummedAmplitudes(k), transverseUpdateFullLabels(k));
            
            end

        end

        %Prunes node with label
        function populationTree = prune(populationTree, label)

            populationTree.root = populationTree.root.prune(label);

        end

        %Update node with label
        function populationTree = updateAmplitudeLabel(populationTree, updateLabel, summedAmplitudes, newLabel)

            populationTree.root = populationTree.root.updateAmplitudeLabel(updateLabel, summedAmplitudes, newLabel);

        end

        function plotTree(populationTree, plotDigits)

            populationTree.root.plotNode(plotDigits);
        
        end

    end    

end