classdef populationTree

    properties (Access = public)

        root longitudinalPopulationNode;
        height double;
        a double;
        TR double;
        f double;
        f_eval double;
        n_tot double;
        hyperpolarizationFactor double;
        pathwayLabelFontSize double;
        amplitudeLabelFontSize double;
        labelOverlapThreshold double;
        n_steady_state double;
        summedTransverseAmplitudes sym;

    end
        
    methods
        
        %Constructor
        function populationTree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, pathwayLabelFontSize, amplitudeLabelFontSize, labelOverlapThreshold, n_steady_state)

            if nargin > 1
            
                populationTree.root = root;
                populationTree.a = a;
                populationTree.TR = TR;
                populationTree.f = f;
                populationTree.f_eval = f_eval;
                populationTree.hyperpolarizationFactor = hyperpolarizationFactor;
                populationTree.n_tot = n_tot+1; %counting a/2 preparation pulse       
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
                populationTree.pathwayLabelFontSize = 12;               
                populationTree.amplitudeLabelFontSize = 12;
                populationTree.labelOverlapThreshold = 0.1;
                populationTree.n_steady_state = inf;

            end

            populationTree.height = 0;
            populationTree.summedTransverseAmplitudes = sym(zeros(1, populationTree.n_tot)); %sum of transverse amplitudes after k pulses 
    
        end

        %Applies pulse to population tree
        function [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject] = applyPulses(populationTreeObject)

            syms x T2p T2 w real;

            for k = 1:min(populationTreeObject.n_tot, populationTreeObject.n_steady_state)

                disp(" ");
     
                if populationTreeObject.height == 0
    
                    [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a/2, populationTreeObject.TR, populationTreeObject.f, populationTreeObject.height); 
    
                elseif populationTreeObject.height == populationTreeObject.n_tot-1
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, populationTreeObject.f_eval, populationTreeObject.height);                
               
                else
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, 1, populationTreeObject.height);  
    
                end
    
                populationTreeObject.height = populationTreeObject.height+1;

                %populationTreeObject = populationTreeObject.pruneMerge(transverseBottomNodes, longitudinalBottomNodes);

                %Calculate sum of transverse amplitudes
                summedAmplitudes = sym(0);
                
                for m = 1:length(transverseBottomNodes)

                    node = transverseBottomNodes(m);
                    
                    partialSpinEchoLocation = -node.coherenceDegreeDirectlyAfterPulse/(360*w); %measured on scale with t0=0 directly after last pulse

                    T2pRelaxation = exp(-abs(x-partialSpinEchoLocation)/T2p);
                    T2Relaxation = exp(-x/T2);
                    phase = exp(1i*360*w*x);
                    summedAmplitudes = summedAmplitudes+node.amplitudeDirectlyAfterPulseWithoutT2p*T2pRelaxation*T2Relaxation*phase; 
                    
                    if k ~= 1
                        disp(newline+"Pathway "+node.label+" has...");
                        disp("  Amplitude (directly after pulse, with T1 and T2 relaxation at t0 = 0, but without any prior T2p relaxation): "+string(node.amplitudeDirectlyAfterPulseWithoutT2p));
                        disp("  T2' Relaxation: "+string(T2pRelaxation));
                        disp("  Additional T2 Relaxation after t = t0 + x: "+string(T2Relaxation));
                        disp("  Phase: "+string(phase));
                        disp("1 of the "+num2str(length(transverseBottomNodes))+" summands for the fitting model for "+num2str(k)+" pulses is the product of the above terms.");
                        disp("Parameters: alpha = "+num2str(populationTreeObject.a)+" deg, TR = "+num2str(populationTreeObject.TR)+" ms, f = "+num2str(populationTreeObject.f));
                    end

                end
  
                populationTreeObject.summedTransverseAmplitudes(k) = summedAmplitudes;

            end

            populationTreeObject.summedTransverseAmplitudes = populationTreeObject.summedTransverseAmplitudes(2:end);

        end

        %Sums branches with the same coherenceDegree by assigning the
        %summed amplitudes to one of the branches, pruning the rest of
        %the branches and updating labels
        function populationTreeObject = pruneMerge(populationTreeObject, transverseBottomNodes, longitudinalBottomNodes)

            disp(" ");

            transverseLabels = strrep(string(zeros(length(transverseBottomNodes), 1)), "0", "");
            transverseAmplitudes = sym(zeros(length(transverseBottomNodes), 1));
            transverseAmplitudeLabels = sym(zeros(length(transverseBottomNodes), 1));
            transversecoherenceDegrees = zeros(length(transverseBottomNodes), 1);

            for k = 1:length(transverseBottomNodes)

                node = transverseBottomNodes(k);

                transverseLabels(k,1) = node.label;
                transverseAmplitudes(k,1) = node.amplitude;
                transverseAmplitudeLabels(k,1) = node.amplitudeLabel;
                transversecoherenceDegrees(k,1) = node.coherenceDegreeDirectlyAfterPulse;

            end

            longitudinalLabels = strrep(string(zeros(length(longitudinalBottomNodes), 1)), "0", "");
            longitudinalAmplitudes = sym(zeros(length(longitudinalBottomNodes), 1));
            longitudinalAmplitudeLabels = sym(zeros(length(longitudinalBottomNodes), 1));
            longitudinalcoherenceDegrees = zeros(length(longitudinalBottomNodes), 1);

            for k = 1:length(longitudinalBottomNodes)

                node = longitudinalBottomNodes(k);

                longitudinalLabels(k,1) = node.label;
                longitudinalAmplitudes(k,1) = node.amplitude;
                longitudinalAmplitudeLabels(k,1) = node.amplitudeLabel;
                longitudinalcoherenceDegrees(k,1) = node.coherenceDegreeDirectlyAfterPulse;

            end
            
            [longitudinalUpdateIndices, longitudinalPruneIndices, longitudinalSummedAmplitudes, longitudinalSummedAmplitudeLabels, longitudinalUpdateFullLabels] = populationTreeMergePruneHelper(longitudinalcoherenceDegrees, longitudinalAmplitudes, longitudinalAmplitudeLabels, longitudinalLabels);
            longitudinalPruneLabels = longitudinalLabels(longitudinalPruneIndices);
            longitudinalUpdateLabels = longitudinalLabels(longitudinalUpdateIndices);

            [transverseUpdateIndices, transversePruneIndices, sTransverseAmplitudes, transverseSummedAmplitudeLabels, transverseUpdateFullLabels] = populationTreeMergePruneHelper(transversecoherenceDegrees, transverseAmplitudes, transverseAmplitudeLabels, transverseLabels);
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
            title("Spin pathways for "+num2str(populationTreeObject.n_tot-1)+" pulses with initial $\frac{\alpha}{2}$ pulse spacing "+num2str(populationTreeObject.f)+" $T_R$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
            xlabel("$t$ (ms)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
            ylabel("Coherence Degree (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

            hold on;
            populationTreeObject.root.plotPathway(populationTreeObject.TR, populationTreeObject.f);

            hold on;
            populationTreeObject.root.plotNode(populationTreeObject.pathwayLabelFontSize, populationTreeObject.amplitudeLabelFontSize, [], [], populationTreeObject.height, populationTreeObject.f, populationTreeObject.labelOverlapThreshold);

            hold on;

            h(2) = plot(0, 0, '.', 'color', [0.6350 0.0780 0.1840], 'MarkerSize', 15, 'DisplayName', '\color[rgb]{0.6350 0.0780 0.1840} Transverse population');
            h(3) = plot(0, 0, '.', 'color', [0.8196 0 1.0000], 'MarkerSize', 15, 'DisplayName', '\color[rgb]{0.8196 0 1.0000} Overlapping transverse & longitudinal populations'); 
            h(4) = line([0, 0], [0, 0], 'color', [0 0.4470 0.7410], 'LineStyle', ':', 'LineWidth', 0.5, 'DisplayName', 'Longitudinal pathway');
            h(5) = line([0, 0], [0, 0], 'color', [0.6350 0.0780 0.1840], 'LineStyle', '--', 'LineWidth', 0.5, 'DisplayName', 'Transverse pathway');
            h(1) = plot(0, 0, '.', 'color', [0 0.4470 0.7410], 'MarkerSize', 15, 'DisplayName', '\color[rgb]{0, 0.4470, 0.7410} Longitudinal population');

            legend(h, 'Location', 'northwest');

            saveas(fig, pwd+"\spinPathways"+num2str(populationTreeObject.n_tot-1)+".fig");
            saveas(fig, pwd+"\spinPathways"+num2str(populationTreeObject.n_tot-1)+".svg");

            close(fig);
        
        end

    end    

end