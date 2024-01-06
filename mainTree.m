close all;
clear all;
clc;

plotDigits=4;

n_tot=3;
Meq=1;
a=50; %deg
TR=10; %ms
f=1/3;

tree=calculatePopulations(n_tot,Meq,a,TR,f);
tree.plotTree(plotDigits);