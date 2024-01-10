close all;
clear all;
clc;
%TODO: lbel with , statt _

%-------------SETTINGS-------------

%Number of pulses
n_tot = 4;

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

%Amplitude label offset (only affects plot)
textOffsetLabelX = -1;
textOffsetLabelY = 0.5;

%Pathway label offset (only affects plot)
textOffsetPathwayX = 0.25;
textOffsetPathwayY = 0.5;

%Pathway label fontsize (only affects plot)
pathwayLabelFontSize = 12;

%Amplitude label fontsize (only affects plot)
amplitudeLabelFontSize = 8;

%-------------END OF SETTINGS-------------

%Create tree with equilibrium magnetization as root
syms M_eq;
root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, 0, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq);
tree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, yScale, textOffsetLabelX, textOffsetLabelY, textOffsetPathwayX, textOffsetPathwayY, pathwayLabelFontSize, amplitudeLabelFontSize);

%Apply pulses
[transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulses();

%Plot tree
tree.plotTree();