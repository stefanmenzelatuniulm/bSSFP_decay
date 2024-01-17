%Sums the transverse amplitudes directly after n_tot pulses

%a: alternating flip angle
%f: initial a/2 pulse spacing
%n_steady_state: calculation stops after n_steady_state pulses because
%steady-state is assumed
%hyperpolarizationFactor: Meq is initially higher by this factor
%createPlot: plot tree?

function summedTransverseAmplitudes = sumTransverseAmplitudes(n_tot, a, TR, f, n_steady_state, hyperpolarizationFactor, createPlot, w0)
    
    %Pathway label fontsize (only affects plot)
    pathwayLabelFontSize = 4;
    
    %Amplitude label fontsize (only affects plot)
    amplitudeLabelFontSize = 5;
    
    %Label overlap threshold (longitudinal amplitude labels are put below
    %pathways if abs(f-0.5)<labelOverlapThreshold) (only affects plot)
    labelOverlapThreshold = 0.1;
    
    %Evaluate tree at f_eval*TR after the last pulse
    f_eval = 1;
    
    %Create tree with equilibrium magnetization as root
    syms M_eq real;
    %parent, transverseChild, longitudinalChild, label, totalTime, coherenceDegree, amplitude, amplitudeLabel, amplitudeDirectlyAfterPulse, amplitudeWithoutT2p, amplitudeDirectlyAfterPulseWithoutT2p, coherenceDegreeDirectlyAfterPulse
    root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq, 0, 0, 0);
    tree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, pathwayLabelFontSize, amplitudeLabelFontSize, labelOverlapThreshold, n_steady_state, w0);
    
    %Apply pulses
    [~, ~, tree] = tree.applyPulses();
    
    %Plot tree
    if createPlot
        tree.plotTree();
    end
    
    summedTransverseAmplitudes = tree.summedTransverseAmplitudes;

end