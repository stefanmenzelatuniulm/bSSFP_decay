%Calculates populations directly after n_tot> = 0 pulses (not counting a/2
%preparation pulse) with equilibrium magnetization Meq, flip angle a, 
%repetition time TR
function tree = calculatePopulations(a, TR, f, f_start, f_eval_end, n_tot, Meq, yScale)
    
    %Create tree with equilibrium magnetization as root
    root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, Meq, 0, sym(1), sym(1));
    tree = populationTree(root, Meq, yScale);

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