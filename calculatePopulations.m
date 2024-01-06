%Calculates populations directly after n_tot>=0 pulses (not counting a/2
%preparation pulse) with equilibrium magnetization Meq, flip angle a,
%repetition time TR
function tree=calculatePopulations(n_tot,Meq,a,TR,f)
    
    %Create tree with equilibrium magnetization as root
    p=population(0,Meq);
    root=transversePopulationNode(p);
    tree=populationTree(root);

    for k=1:n_tot
    
        [transverseBottomNodes,longitudinalBottomNodes,tree]=tree.applyPulse(a,TR,f);
        tree=tree.pruneMerge(transverseBottomNodes,longitudinalBottomNodes);
    
    end

end