close all;
clear all;
clc;

%TODO: return value for prune, for feedback -> verbosity
%Better comments and formatting
%Sum and return bottom transverse amplitudes

%-------------SETTINGS-------------

%Number of pulses
n_tot = 3;

%Flip angle (deg)
a = 50;

%Repetition time (ms)
TR = 10;

%Initial alpha/2 pulse spacing f*TR
f = 1/2;

%Evaluation at f_eval*TR after the last pulse
f_eval = 1;

%Hyperpolarization factor
hyperpolarizationFactor = 1;

%y axis scale (only affects plot)
yScale = 1;

%Pathway label fontsize (only affects plot)
pathwayLabelFontSize = 14;

%Amplitude label fontsize (only affects plot)
amplitudeLabelFontSize = 10;

%Label overlap threshold (longitudinal amplitude labels are put below
%pathways if abs(f-0.5)<labelOverlapThreshold)
labelOverlapThreshold = 0.1;

%-------------END OF SETTINGS-------------

%Create tree with equilibrium magnetization as root
syms M_eq;
root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, 0, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq);
tree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, yScale, pathwayLabelFontSize, amplitudeLabelFontSize, labelOverlapThreshold);

%Apply pulses
[transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulses();

%Plot tree
tree.plotTree();