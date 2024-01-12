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
        pathwayLabelFontSize double;
        amplitudeLabelFontSize double;
        labelOverlapThreshold double;
        n_steady_state int64;
        summedTransverseAmplitudes sym;

    end
        
    methods
        
        %Constructor
        function populationTree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, yScale, pathwayLabelFontSize, amplitudeLabelFontSize, labelOverlapThreshold, n_steady_state)

            if nargin > 1
            
                populationTree.root = root;
                populationTree.a = a;
                populationTree.TR = TR;
                populationTree.f = f;
                populationTree.f_eval = f_eval;
                populationTree.hyperpolarizationFactor = hyperpolarizationFactor;
                populationTree.n_tot = n_tot+1; %counting a/2 preparation pulse
                populationTree.yScale = yScale;             
                populationTree.pathwayLabelFontSize = pathwayLabelFontSize;
                populationTree.amplitudeLabelFontSize = amplitudeLabelFontSize;
                populationTree.labelOverlapThreshold = labelOverlapThreshold;
                populationTree.n_steady_state = n_steady_state;

            else

                populationTree.root = emptyNode();
                populationTree.a = 90;
                populationTree.TR = 10;
                populationTree.f = 0.5;
                populationTree.f_eval = 0.5;
                populationTree.n_tot = 1;
                populationTree.hyperpolarizationFactor = 1;
                populationTree.yScale = 1;
                populationTree.pathwayLabelFontSize = 12;               
                populationTree.amplitudeLabelFontSize = 12;
                populationTree.labelOverlapThreshold = 0.1;
                populationTree.n_steady_state = intmax("int64");

            end

            populationTree.height = 0;
            populationTree.summedTransverseAmplitudes = sym(zeros(1, n_tot)); %sum of transverse amplitudes after k pulses 
    
        end

        %Applies pulse to population tree
        function [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject] = applyPulses(populationTreeObject)

            syms x T2p T2;

            for k = 1:min(populationTreeObject.n_tot, populationTreeObject.n_steady_state)
     
                if populationTreeObject.height == 0
    
                    [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a/2, populationTreeObject.TR, populationTreeObject.f, populationTreeObject.height, populationTreeObject.yScale); %bottomNodes is list with entries: label+"#"+string(amplitude)+"#"+string(amplitudeLabel)+"#"+string(dephasingDegree)
    
                elseif populationTreeObject.height == populationTreeObject.n_tot-1
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, populationTreeObject.f_eval, populationTreeObject.height, populationTreeObject.yScale);                
               
                else
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, 1, populationTreeObject.height, populationTreeObject.yScale);  
    
                end
    
                populationTreeObject.height = populationTreeObject.height+1;

                %Calculate sum of transverse amplitudes
                summedAmplitudes = sym(0);
                
                for m = 1:length(transverseBottomNodes)
                s
                    node = transverseBottomNodes(m);
                    
                    summedAmplitudes = summedAmplitudes+node.amplitudeDirectlyAfterPulse*exp(-abs(x-abs(node.dephasingDegree))/T2p)*exp(abs(node.dephasingDegree)/T2p)*exp(-x/T2); 
                
                end
    
                summedAmplitudes = simplify(summedAmplitudes, "IgnoreAnalyticConstraints", true);

                populationTreeObject.summedTransverseAmplitudes(k) = summedAmplitudes;
               
                populationTreeObject = populationTreeObject.pruneMerge(transverseBottomNodes, longitudinalBottomNodes);

            end

            populationTreeObject.summedTransverseAmplitudes = populationTreeObject.summedTransverseAmplitudes(2:end);

        end

        %Sums branches with the same dephasing degree by assigning the
        %summed amplitudes to one of the branches, pruning the rest of
        %the branches and updating labels
        function populationTreeObject = pruneMerge(populationTreeObject, transverseBottomNodes, longitudinalBottomNodes)

            transverseLabels = strrep(string(zeros(length(transverseBottomNodes), 1)), "0", "");
            transverseAmplitudes = sym(zeros(length(transverseBottomNodes), 1));
            transverseAmplitudeLabels = sym(zeros(length(transverseBottomNodes), 1));
            transverseDephasingDegrees = zeros(length(transverseBottomNodes), 1);

            for k = 1:length(transverseBottomNodes)

                node = transverseBottomNodes(k);

                transverseLabels(k,1) = node.label;
                transverseAmplitudes(k,1) = node.amplitude;
                transverseAmplitudeLabels(k,1) = node.amplitudeLabel;
                transverseDephasingDegrees(k,1) = node.dephasingDegree;

            end

            longitudinalLabels = strrep(string(zeros(length(longitudinalBottomNodes), 1)), "0", "");
            longitudinalAmplitudes = sym(zeros(length(longitudinalBottomNodes), 1));
            longitudinalAmplitudeLabels = sym(zeros(length(longitudinalBottomNodes), 1));
            longitudinalDephasingDegrees = zeros(length(longitudinalBottomNodes), 1);

            for k = 1:length(longitudinalBottomNodes)

                node = longitudinalBottomNodes(k);

                longitudinalLabels(k,1) = node.label;
                longitudinalAmplitudes(k,1) = node.amplitude;
                longitudinalAmplitudeLabels(k,1) = node.amplitudeLabel;
                longitudinalDephasingDegrees(k,1) = node.dephasingDegree;

            end
            
            [longitudinalUpdateIndices, longitudinalPruneIndices, longitudinalSummedAmplitudes, longitudinalSummedAmplitudeLabels, longitudinalUpdateFullLabels] = populationTreeMergePruneHelper(longitudinalDephasingDegrees, longitudinalAmplitudes, longitudinalAmplitudeLabels, longitudinalLabels);
            longitudinalPruneLabels = longitudinalLabels(longitudinalPruneIndices);
            longitudinalUpdateLabels = longitudinalLabels(longitudinalUpdateIndices);

            [transverseUpdateIndices, transversePruneIndices, sTransverseAmplitudes, transverseSummedAmplitudeLabels, transverseUpdateFullLabels] = populationTreeMergePruneHelper(transverseDephasingDegrees, transverseAmplitudes, transverseAmplitudeLabels, transverseLabels);
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
            
                populationTreeObject = populationTreeObject.updateAmplitudeLabel(transverseUpdateLabels(k), sTransverseAmplitudes(k), transverseSummedAmplitudeLabels(k), transverseUpdateFullLabels(k));
            
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

            fig = figure('WindowState','maximized');
            ax = gca;
            ax.FontSize = 14;
            title("Spin pathways for "+num2str(populationTreeObject.n_tot)+" pulses with initial $\frac{\alpha}{2}$ pulse spacing "+num2str(populationTreeObject.f)+" $T_R$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
            xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
            ylabel("Dephasing degree (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

            hold on;
            populationTreeObject.root.plotPathway(populationTreeObject.TR, populationTreeObject.f);

            hold on;
            populationTreeObject.root.plotNode(populationTreeObject.pathwayLabelFontSize, populationTreeObject.amplitudeLabelFontSize, [], [], populationTreeObject.height, populationTreeObject.f, populationTreeObject.labelOverlapThreshold);

            hold on;

            h(2) = plot(0, 0, '.', 'color', [0.6350 0.0780 0.1840], 'MarkerSize', 20, 'DisplayName', '\color[rgb]{0.6350 0.0780 0.1840} Transverse population');
            h(3) = plot(0, 0, '.', 'color', [0.8196 0 1.0000], 'MarkerSize', 20, 'DisplayName', '\color[rgb]{0.8196 0 1.0000} Overlapping transverse & longitudinal populations'); 
            h(4) = line([0, 0], [0, 0], 'color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 0.5, 'DisplayName', 'Longitudinal pathway');
            h(5) = line([0, 0], [0, 0], 'color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 0.5, 'DisplayName', 'Transverse pathway');
            h(1) = plot(0, 0, '.', 'color', [0 0.4470 0.7410], 'MarkerSize', 20, 'DisplayName', '\color[rgb]{0, 0.4470, 0.7410} Longitudinal population');

            legend(h, 'Location', 'northwest');

            yExtentMax = populationTreeObject.TR*populationTreeObject.yScale*populationTreeObject.f+populationTreeObject.TR*populationTreeObject.yScale*(double(populationTreeObject.height)-2)+populationTreeObject.f_eval*populationTreeObject.TR*populationTreeObject.yScale*(double(populationTreeObject.height)-2);
            yExtentMin = max(0, populationTreeObject.TR*populationTreeObject.yScale*populationTreeObject.f+populationTreeObject.TR*populationTreeObject.yScale*(double(populationTreeObject.height)-2));
            ylim([-max(0.2, yExtentMin), yExtentMax]*1.15);
            xl = xlim;
            xlim([xl(1), xl(2)+0.25*(double(populationTreeObject.height-2)*populationTreeObject.TR+populationTreeObject.f_eval*populationTreeObject.TR)]);

            saveas(fig, pwd+"\spinPathways"+num2str(populationTreeObject.n_tot-1)+".fig");
            saveas(fig, pwd+"\spinPathways"+num2str(populationTreeObject.n_tot-1)+".svg");

            close(fig);
        
        end

    end    

end