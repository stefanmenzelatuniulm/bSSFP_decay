%Calculates populations directly after n_tot> = 0 pulses (not counting a/2
%preparation pulse) with equilibrium magnetization Meq, flip angle a, 
%repetition time TR
function tree = calculatePopulations(a, TR, f, f_start, f_eval_end, n_tot, yScale, hyperpolarizationFactor)
    
    %Create tree with equilibrium magnetization as root
    syms Meq;
    root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, 0, hyperpolarizationFactor*Meq, hyperpolarizationFactor*Meq);
    tree = populationTree(root, yScale);

    for k = 1:n_tot
    
        if k == 1
            [transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulse(a, TR, f_start);
        elseif k == n_tot
            [transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulse(a, TR, f_eval_end);
        else
            [transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulse(a, TR, f);
        end

        tree = tree.pruneMerge(transverseBottomNodes, longitudinalBottomNodes);
    
    end

end