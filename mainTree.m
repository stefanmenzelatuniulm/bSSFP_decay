close all;
clear all;
clc;

%-------------SETTINGS-------------

%Number of pulses
n_tot = 1;

%Flip angle (deg)
a = 50;

%Repetition time (ms)
TR = 10;

%Initial alpha/2 pulse spacing f*TR
f = 1/3;

%Evaluation at f_eval*TR after the last pulse
f_eval = 1;

%Hyperpolarization factor
hyperpolarizationFactor = 1;

%y axis scale
yScale = 1;

%-------------END OF SETTINGS-------------

%Create tree with equilibrium magnetization as root
syms M_eq;
root = longitudinalPopulationNode(emptyNode(), emptyNode(), emptyNode(), "", 0, 0, 0, hyperpolarizationFactor*M_eq, hyperpolarizationFactor*M_eq);
tree = populationTree(root, a, TR, f, f_eval, n_tot, hyperpolarizationFactor, yScale);

%Apply pulses
[transverseBottomNodes, longitudinalBottomNodes, tree] = tree.applyPulses();

%Plot tree
tree.plotTree();