close all;
clear all;
clc;

plotDigits = 4;

a = 50; %deg
TR = 10; %ms
f = 1;
f_start = 1/3;
f_eval_end = 1;
n_tot = 3;
Meq = 1;

yScale = 1;

tree = calculatePopulations(a, TR, f, f_start, f_eval_end, n_tot, Meq, yScale);
tree.plotTree(plotDigits);