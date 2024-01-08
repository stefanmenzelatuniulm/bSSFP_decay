classdef populationTree

    properties (Access = public)

        root longitudinalPopulationNode;
        height int64;
        a double;
        TR double;
        f double;
        f_eval double;
        n_tot int64;
        hyperpolarizationFactor double;
        yScale double;

    end
        
    methods
        
        %Constructor
        function populationTree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, yScale)

            if nargin > 1
            
                populationTree.root = root;
                populationTree.a = a;
                populationTree.TR = TR;
                populationTree.f = f;
                populationTree.f_eval = f_eval;
                populationTree.hyperpolarizationFactor = hyperpolarizationFactor;
                populationTree.n_tot = n_tot;
                populationTree.yScale = yScale;

            else

                populationTree.root = emptyNode();
                populationTree.a = 90;
                populationTree.TR = 10;
                populationTree.f = 0.5;
                populationTree.f_eval = 0.5;
                populationTree.n_tot = 1;
                populationTree.hyperpolarizationFactor = 1;
                populationTree.yScale = 1;

            end

            populationTree.height = 0;
    
        end

        %Applies pulse to population tree
        function [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject] = applyPulses(populationTreeObject)

            for k = 1:populationTreeObject.n_tot
            
                if populationTreeObject.height == 0
    
                    [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a, populationTreeObject.TR, populationTreeObject.f, populationTreeObject.height, populationTreeObject.yScale); %bottomNodes is list with entries: label+"#"+string(amplitude)+"#"+string(amplitudeLabel)+"#"+string(dephasingDegree)
    
                elseif populationTreeObject.height == populationTreeObject.n_tot
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a, populationTreeObject.TR, populationTreeObject.f_eval, populationTreeObject.height, populationTreeObject.yScale);
    
                else
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a, populationTreeObject.TR, 1, populationTreeObject.height, populationTreeObject.yScale);  
    
                end
    
                populationTreeObject.height = populationTreeObject.height+1;
                populationTreeObject = populationTreeObject.pruneMerge(transverseBottomNodes, longitudinalBottomNodes);

            end

        end

        %Sums branches with the same dephasing degree by assigning the
        %summed amplitudes to one of the branches, pruning the rest of
        %the branches and updating labels
        function populationTreeObject = pruneMerge(populationTreeObject, transverseBottomNodes, longitudinalBottomNodes)

            transverseBottomNodes = split(transpose(transverseBottomNodes), "#");
            longitudinalBottomNodes = split(transpose(longitudinalBottomNodes), "#");

            dim = size(transverseBottomNodes);

            if dim(2) == 1
                transverseBottomNodes = transpose(transverseBottomNodes);
                longitudinalBottomNodes = transpose(longitudinalBottomNodes);
            end
            
            transverseLabels = transverseBottomNodes(:, 1);
            transverseAmplitudes = str2sym(transverseBottomNodes(:, 2));
            transverseAmplitudeLabels = str2sym(transverseBottomNodes(:, 3));
            transverseDephasingDegrees = double(transverseBottomNodes(:, 4));

            longitudinalLabels = longitudinalBottomNodes(:, 1);
            longitudinalAmplitudes = str2sym(longitudinalBottomNodes(:, 2));
            longitudinalAmplitudeLabels = str2sym(longitudinalBottomNodes(:, 3));
            longitudinalDephasingDegrees = double(longitudinalBottomNodes(:, 4));
            
            [longitudinalUpdateIndices, longitudinalPruneIndices, longitudinalSummedAmplitudes, longitudinalSummedAmplitudeLabels, longitudinalUpdateFullLabels] = populationTreeMergePruneHelper(longitudinalDephasingDegrees, longitudinalAmplitudes, longitudinalAmplitudeLabels, longitudinalLabels);
            longitudinalPruneLabels = longitudinalLabels(longitudinalPruneIndices);
            longitudinalUpdateLabels = longitudinalLabels(longitudinalUpdateIndices);

            [transverseUpdateIndices, transversePruneIndices, transverseSummedAmplitudes, transverseSummedAmplitudeLabels, transverseUpdateFullLabels] = populationTreeMergePruneHelper(transverseDephasingDegrees, transverseAmplitudes, transverseAmplitudeLabels, transverseLabels);
            transversePruneLabels = transverseLabels(transversePruneIndices);
            transverseUpdateLabels = transverseLabels(transverseUpdateIndices);

            for k = 1:length(longitudinalPruneLabels)
            
                populationTreeObject = populationTreeObject.prune(longitudinalPruneLabels(k));
            
            end

            for k = 1:length(longitudinalUpdateLabels)
                populationTreeObject = populationTreeObject.updateAmplitudeLabel(longitudinalUpdateLabels(k), longitudinalSummedAmplitudes(k), longitudinalSummedAmplitudeLabels(k), longitudinalUpdateFullLabels(k));
            
            end

            for k = 1:length(transversePruneLabels)
            
                populationTreeObject = populationTreeObject.prune(transversePruneLabels(k));
            
            end

            for k = 1:length(transverseUpdateLabels)
            
                populationTreeObject = populationTreeObject.updateAmplitudeLabel(transverseUpdateLabels(k), transverseSummedAmplitudes(k), transverseSummedAmplitudeLabels(k), transverseUpdateFullLabels(k));
            
            end

        end

        %Prunes node with label
        function populationTreeObject = prune(populationTreeObject, label)

            populationTreeObject.root = populationTreeObject.root.prune(label);

        end

        %Update node with label
        function populationTreeObject = updateAmplitudeLabel(populationTreeObject, updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel)

            populationTreeObject.root = populationTreeObject.root.updateAmplitudeLabel(updateLabel, summedAmplitudes, summedAmplitudeLabels, newLabel);

        end

        %Plot tree
        function plotTree(populationTreeObject)

            fig=figure('WindowState','maximized');
            ax = gca;
            ax.FontSize = 14;
            title("Spin pathways for "+num2str(populationTreeObject.n_tot)+" pulses with $\alpha = $ "+num2str(populationTreeObject.a)+" $^{\circ}$, $T_R = $ "+num2str(populationTreeObject.TR)+" ms, initial $\frac{\alpha}{2}$ pulse spacing "+num2str(populationTreeObject.f)+" $T_R$, hyperpolarization factor "+num2str(populationTreeObject.hyperpolarizationFactor)+", evaluated at "+num2str(populationTreeObject.f_eval)+" $T_R$","interpreter","latex",'fontweight','bold','fontsize',14);
            xlabel("$t$ (ms)","interpreter","latex",'fontweight','bold','fontsize',14);
            ylabel("Dephasing degree (ms)","interpreter","latex",'fontweight','bold','fontsize',14);

            hold on;
            populationTreeObject.root.plotNode();

            hold on;

            h(6) = line([0, 0], [0, 0], 'color', [0.2 0.2 0.2], 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', 'Longitudinal pathway');
            h(5) = line([0, 0], [0, 0], 'color', [0.2 0.2 0.2], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Transverse pathway');
            h(3) = plot(0, 0, '.', 'Color', [0.4940 0.1840 0.5560], 'MarkerSize', 30, 'DisplayName', '\color[rgb]{0.4940, 0.1840, 0.5560} Pathway label');
            h(4) = plot(0, 0, '.', 'Color', [0.8500 0.3250 0.0980], 'MarkerSize', 30, 'DisplayName', '\color[rgb]{0.8500, 0.3250, 0.0980} Pathway amplitude');  
            h(2) = plot(0, 0, '.', 'color', [0.4660 0.6740 0.1880], 'MarkerSize', 30, 'DisplayName', '\color[rgb]{0.4660, 0.6740, 0.1880} Transverse population');
            h(1) = plot(0, 0, '.', 'color', [0 0.4470 0.7410], 'MarkerSize', 30, 'DisplayName', '\color[rgb]{0, 0.4470, 0.7410} Longitudinal population');

            legend(h, 'Location', 'northwest');

            saveas(fig,pwd+"\spinPathways.fig");
            saveas(fig,pwd+"\spinPathways.png");

            close(fig);
        
        end

    end    

end