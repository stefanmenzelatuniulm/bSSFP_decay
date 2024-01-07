close all;
clear all;
clc;

a = 50;
TR = 10;
f = 1;
f_start = 1/3;
f_eval_end = 1;
n_tot = 3;
hyperpolarizationFactor = 1;

yScale = 1;
plotDigits = 4;

tree = calculatePopulations(a, TR, f, f_start, f_eval_end, n_tot, yScale, hyperpolarizationFactor);
tree.plotTree(plotDigits);