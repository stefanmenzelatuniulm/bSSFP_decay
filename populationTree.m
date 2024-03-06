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
        n_steady_state double;
        summedTransverseAmplitudes sym;
        w0 double;
        T2sDominated logical;
        decreaseVerbosity logical;

    end
        
    methods
        
        %Constructor
        function populationTree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, n_steady_state, w0, T2sDominated, decreaseVerbosity)

            if nargin > 1
            
                populationTree.root = root;
                populationTree.a = a;
                populationTree.TR = TR;
                populationTree.f = f;
                populationTree.f_eval = f_eval;
                populationTree.hyperpolarizationFactor = hyperpolarizationFactor;
                populationTree.n_tot = n_tot+1; %counting a/2 preparation pulse       
                populationTree.n_steady_state = n_steady_state+1; %counting a/2 preparation pulse
                populationTree.w0 = w0;
                populationTree.T2sDominated = T2sDominated;
                populationTree.decreaseVerbosity = decreaseVerbosity;

            else

                populationTree.root = emptyNode();
                populationTree.a = 90;
                populationTree.TR = 10;
                populationTree.f = 0.5;
                populationTree.f_eval = 0.5;
                populationTree.n_tot = 1;
                populationTree.hyperpolarizationFactor = 1;
                populationTree.n_steady_state = inf;
                populationTree.w0 = 0;
                populationTree.T2sDominated = false;
                populationTree.decreaseVerbosity = false;

            end

            populationTree.height = 0;
            populationTree.summedTransverseAmplitudes = sym(zeros(1, populationTree.n_tot)); %sum of transverse amplitudes after k pulses 
    
        end

        %Applies pulse to population tree
        function [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject] = applyPulses(populationTreeObject)
            
            syms x Psi real;
            T2s = sym("T2s", "real");
            if ~populationTreeObject.T2sDominated
                syms T2;
            end

            for k = 1:min(populationTreeObject.n_tot, populationTreeObject.n_steady_state)

                if k>1
                    disp(" ");
                    disp("Applying pulse "+num2str(k-1)+"/"+num2str(populationTreeObject.n_tot-1));
                else
                    disp("Applying initial alpha/2 pulse (counted as 0th pulse)");
                end
     
                if populationTreeObject.height ==  0
    
                    [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(populationTreeObject.a/2, populationTreeObject.TR, populationTreeObject.f, populationTreeObject.height, populationTreeObject.T2sDominated, populationTreeObject.decreaseVerbosity); 
    
                elseif populationTreeObject.height ==  populationTreeObject.n_tot-1
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, populationTreeObject.f_eval, populationTreeObject.height, populationTreeObject.T2sDominated, populationTreeObject.decreaseVerbosity);                
               
                else
    
                     [transverseBottomNodes, longitudinalBottomNodes, populationTreeObject.root] = populationTreeObject.root.applyPulse(((-1)^(k+1))*populationTreeObject.a, populationTreeObject.TR, 1, populationTreeObject.height, populationTreeObject.T2sDominated, populationTreeObject.decreaseVerbosity);  
    
                end
    
                populationTreeObject.height = populationTreeObject.height+1;

                if k~=1
                    disp(" ");
                    disp("Preparing to prune/merge overlapping pathways from pulse "+num2str(k-1)+"/"+num2str(populationTreeObject.n_tot-1));
                end

                populationTreeObject = populationTreeObject.pruneMerge(transverseBottomNodes, longitudinalBottomNodes);

                %Calculate sum of transverse amplitudes
                if k~=1
                    disp("Summing "+num2str(length(transverseBottomNodes))+" transverse pathways from pulse "+num2str(k-1)+"/"+num2str(populationTreeObject.n_tot-1));
                end
                summedAmplitudes = sym(0);
                
                for m = 1:length(transverseBottomNodes)

                    node = transverseBottomNodes(m);

                    if ~populationTreeObject.T2sDominated
                        T2Relaxation = exp(-x/T2);
                    else
                        T2Relaxation = 1;
                    end

                    phase = exp(1i*2*pi*(populationTreeObject.w0+Psi)*(x+node.dephasingTimeDirectlyAfterPulse))*exp(-(1/T2s)*abs(x+node.dephasingTimeDirectlyAfterPulse));
                    summedAmplitudes = summedAmplitudes+node.amplitudeDirectlyAfterPulseWithoutT2s*T2Relaxation*phase; 
                    
                    if k ~=  1

                        if ~populationTreeObject.decreaseVerbosity
                            disp(newline+"  Pathway "+node.label+" has:");
                            if ~populationTreeObject.T2sDominated
                                disp("      Amplitude directly after pulse, with prior T1 and T2 Relaxation up to t0 = 0, but without any prior T2* Relaxation): "+string(node.amplitudeDirectlyAfterPulseWithoutT2s));
                                disp("      Additional T2 Relaxation after t = t0 + x: "+string(T2Relaxation));
                            end
                            disp("      Phase & T2* Relaxation: "+string(phase));
                            disp("  1 of the "+num2str(length(transverseBottomNodes))+" summands for the fitting model for "+num2str(k-1)+" pulses is the product of the first three terms.");
                            disp("  Parameters: alpha = "+num2str(populationTreeObject.a)+" deg, TR = "+num2str(populationTreeObject.TR)+" ms, f = "+num2str(populationTreeObject.f));
                        end

                    end

                end
  
                populationTreeObject.summedTransverseAmplitudes(k) = summedAmplitudes;

            end

            if populationTreeObject.n_tot > populationTreeObject.n_steady_state

                disp(" ");
            
                for k = populationTreeObject.n_steady_state+1:populationTreeObject.n_tot

                    disp("Steady-state assumption reached: Signal after pulse "+num2str(k-1)+" is the same as the signal after pulse "+num2str(populationTreeObject.n_steady_state-1));
                
                    populationTreeObject.summedTransverseAmplitudes(k) = populationTreeObject.summedTransverseAmplitudes(k-1);

                end
            
            end

            populationTreeObject.summedTransverseAmplitudes = populationTreeObject.summedTransverseAmplitudes(2:end);

        end

        %Sums pathways with the same coherenceDegree by assigning the
        %summed amplitudes to one of the pathways, pruning the rest of
        %the pathways and updating labels
        function populationTreeObject = pruneMerge(populationTreeObject, transverseBottomNodes, longitudinalBottomNodes)

            if ~populationTreeObject.decreaseVerbosity || populationTreeObject.height>1
                disp(" ");
            end

            transverseLabels = strrep(string(zeros(length(transverseBottomNodes), 1)), "0", "");
            transversecoherenceDegrees = zeros(length(transverseBottomNodes), 1);
            transverseAmplitudesWithoutT2s = sym(zeros(length(transverseBottomNodes), 1));
            transverseAmplitudesDirectlyAfterPulseWithoutT2s = sym(zeros(length(transverseBottomNodes), 1));

            for k = 1:length(transverseBottomNodes)

                node = transverseBottomNodes(k);

                transverseLabels(k, 1) = node.label;
                transversecoherenceDegrees(k, 1) = node.coherenceDegreeDirectlyAfterPulse;
                transverseAmplitudesWithoutT2s(k, 1) = node.amplitudeWithoutT2s;
                transverseAmplitudesDirectlyAfterPulseWithoutT2s(k, 1) = node.amplitudeDirectlyAfterPulseWithoutT2s;

            end

            longitudinalLabels = strrep(string(zeros(length(longitudinalBottomNodes), 1)), "0", "");
            longitudinalcoherenceDegrees = zeros(length(longitudinalBottomNodes), 1);
            longitudinalAmplitudesWithoutT2s = sym(zeros(length(longitudinalBottomNodes), 1));
            longitudinalAmplitudesDirectlyAfterPulseWithoutT2s = sym(zeros(length(longitudinalBottomNodes), 1));

            for k = 1:length(longitudinalBottomNodes)

                node = longitudinalBottomNodes(k);

                longitudinalLabels(k, 1) = node.label;
                longitudinalcoherenceDegrees(k, 1) = node.coherenceDegreeDirectlyAfterPulse;
                longitudinalAmplitudesWithoutT2s(k, 1) = node.amplitudeWithoutT2s;
                longitudinalAmplitudesDirectlyAfterPulseWithoutT2s(k, 1) = node.amplitudeDirectlyAfterPulseWithoutT2s;                

            end
            
            [longitudinalUpdateIndices, longitudinalPruneIndices, longitudinalUpdateFullLabels, longitudinalSummedAmplitudesWithoutT2s, longitudinalSummedAmplitudesDirectlyAfterPulseWithoutT2s] = populationTreeMergePruneHelper(longitudinalcoherenceDegrees, longitudinalLabels, longitudinalAmplitudesWithoutT2s, longitudinalAmplitudesDirectlyAfterPulseWithoutT2s);
            longitudinalPruneLabels = longitudinalLabels(longitudinalPruneIndices);
            longitudinalUpdateLabels = longitudinalLabels(longitudinalUpdateIndices);

            [transverseUpdateIndices, transversePruneIndices, transverseUpdateFullLabels, transverseSummedAmplitudesWithoutT2s, transverseSummedAmplitudesDirectlyAfterPulseWithoutT2s] = populationTreeMergePruneHelper(transversecoherenceDegrees, transverseLabels, transverseAmplitudesWithoutT2s, transverseAmplitudesDirectlyAfterPulseWithoutT2s);
            transversePruneLabels = transverseLabels(transversePruneIndices);
            transverseUpdateLabels = transverseLabels(transverseUpdateIndices);

            if populationTreeObject.height>1
                disp("  Pruning "+num2str(length(longitudinalPruneLabels))+" overlapping longitudinal pathways");
            end

            for k = 1:length(longitudinalPruneLabels)
            
                populationTreeObject = populationTreeObject.prune(longitudinalPruneLabels(k));
            
            end
            
            longIndices = find(longitudinalUpdateLabels~=longitudinalUpdateFullLabels);

            if ~populationTreeObject.decreaseVerbosity
                disp(" ");
            end
            if populationTreeObject.height>1
                disp("  Updating "+num2str(length(longIndices))+" overlapping longitudinal pathways");
            end

            for k = 1:length(longIndices)

                lind = longIndices(k);
                populationTreeObject = populationTreeObject.updateAmplitudeLabel(longitudinalUpdateLabels(lind), longitudinalUpdateFullLabels(lind), longitudinalSummedAmplitudesWithoutT2s(lind), longitudinalSummedAmplitudesDirectlyAfterPulseWithoutT2s(lind));
            
            end

            if ~populationTreeObject.decreaseVerbosity
                disp(" ");
            end
            if populationTreeObject.height>1
                disp("  Pruning "+num2str(length(transversePruneLabels))+" overlapping transverse pathways");
            end

            for k = 1:length(transversePruneLabels)
            
                populationTreeObject = populationTreeObject.prune(transversePruneLabels(k));
            
            end

            transvIndices = find(transverseUpdateLabels~=transverseUpdateFullLabels);

            if ~populationTreeObject.decreaseVerbosity
                disp(" ");
            end      
            if populationTreeObject.height>1
                disp("  Updating "+num2str(length(transvIndices))+" overlapping transverse pathways");
            end

            for k = 1:length(transvIndices)
            
                tind = transvIndices(k);
                populationTreeObject = populationTreeObject.updateAmplitudeLabel(transverseUpdateLabels(tind), transverseUpdateFullLabels(tind), transverseSummedAmplitudesWithoutT2s(tind), transverseSummedAmplitudesDirectlyAfterPulseWithoutT2s(tind));
            
            end

            if populationTreeObject.decreaseVerbosity && populationTreeObject.height>1
                disp(" ");
            end

        end

        %Prune node with label
        function populationTreeObject = prune(populationTreeObject, label)

            populationTreeObject.root = populationTreeObject.root.prune(label, populationTreeObject.decreaseVerbosity);

        end

        %Update node with label
        function populationTreeObject = updateAmplitudeLabel(populationTreeObject, updateLabel, newLabel, summedAmplitudesWithoutT2s, summedAmplitudesDirectlyAfterPulseWithoutT2s)

            populationTreeObject.root = populationTreeObject.root.updateAmplitudeLabel(updateLabel, newLabel, summedAmplitudesWithoutT2s, summedAmplitudesDirectlyAfterPulseWithoutT2s, populationTreeObject.decreaseVerbosity);

        end

    end    

end