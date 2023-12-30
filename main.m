clear all;
close all;
clc; 

%-------------SETTINGS-------------

%Number of isochromats
ns=1024;
            
%Number of -/+ alpha pulses -1 (no sampling after last pulse due to +/-
%alpha/2 tip-back pulse)
n_tot=32;

%Metabolite properties (13C lactate)
w0=300; %in Hz, rotating frame (gamma*B0 is filtered by heterodyne mixing)
FWHM=21; %in Hz %Simulation depends sensitively on FWHM
T1=13.1; %in s
T2=0.6; %in s

%Range of flip angles in degree
%Calculation much faster for small cosd(alpha) because steady-state approximation can be applied faster 
amin=35;
amax=145;
na=23; %resolution: 5 degree

%Range of repetition times in ms
TRmin=5.5;
TRmax=20.5;
nTR=16; %resolution: 1/10 ms

%Range of times TR*f between initial alpha/2 and first -alpha pulse, if TR
%is the time between -/+ alpha pulses
f=[1/3 0.5];

%Calculate signal at time t_eval=f_eval*TR, measured from the end of the
%pulse train. f_eval=0 e.g. calculates the signal directly after the end of
%the pulse train
f_eval_min=0;
f_eval_max=1;
nf_eval=1001; %resolution: 1/1000 TR

%Equilibrium magnetization
Meq=1;

%splitfactor loop iterations are used in vectorizedM -> high splitfactor
%causes less RAM usage, but vectorization is not as efficient
splitfactor=128; 

%Recalculate M, or read existing M from save file M.txt?
recalculateM=false;

%Recalculate C1, C2, C2s or read from save file C.txt?
recalculateC=false;

%Steady-state approximation: Set (E1*cosd(alpha))^k = 0 if abs((E1*cos(alpha))^k) < epsilon
%-> Calculations much faster if cosd(alpha) is small
epsilon=0.01;

%Upper bounds for T1 and T2 (regarding fit, and steady-state approximation)
T1max=100; %in s
T2max=20; %in s

%-------------END OF SETTINGS-------------

%Flip angles in degree
a=linspace(amin,amax,na);

%Repetition times in s
TR=linspace(TRmin,TRmax,nTR);

%Evaluate at f_eval*TR
f_eval=linspace(f_eval_min,f_eval_max,nf_eval);

%Lorentzian distribution with FWHM around w0 for Larmor frequencies in Hz
%(Lorentzian distribution is special case of Student's t distribution with
%1 degree of freedom) Care: sigma is HWHM
pd=makedist('tLocationScale','mu',w0,'sigma',FWHM/2,'nu',1);
w=random(pd,1,ns);

%Calculate M Cant fully exploit the vectorization of vectorizedM due to RAM
%overflow Balance has to be found of CPU time vs RAM usage
if recalculateM

    M=zeros(na,nTR,length(f),nf_eval,n_tot); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
    
    for k=1:nTR
        
        disp("Calculating M for the "+num2str(k)+"th repetition time...");
        M(:,k,:,:,:)=vectorizedM(a,TR(k),w,f,f_eval,n_tot,Meq,T1,T2,splitfactor); 
    
    end

    %Save M and w
    save(pwd+"\"+"M.mat","M","-v7.3");
    save(pwd+"\"+"w.mat","w","-v7.3");

else

    load(pwd+"\"+"M.mat");
    load(pwd+"\"+"w.mat");

end

mkdir("Figures");
deleteFigures("Figures");

%Plot histogram of w
plotHist(w,w0,FWHM);

%Calculate relaxation coefficients
if recalculateC

    [C1 C2 C2s]=decayCoefficients(a,TR,f,n_tot,epsilon,T1max,Meq);

    %Save C1, C2, C2s
    save(pwd+"\"+"C1.mat","C1","-v7.3");
    save(pwd+"\"+"C2.mat","C2","-v7.3");
    save(pwd+"\"+"C2s.mat","C2s","-v7.3");

else

    load(pwd+"\"+"C1.mat");
    load(pwd+"\"+"C2.mat");
    load(pwd+"\"+"C2s.mat");

end

ft=zeros(na,nTR,length(f),n_tot);

for k=1:length(f)

    for m=1:na

        for o=1:nTR

            for l=2:n_tot %1st pulse: signal does not depend on T1 -> problems fitting -> start at 2nd pulse

                M_=permute(M(m,o,k,:,l),[4 1 2 3 5]);

                ft_=plotM2Dfit(M_,f_eval,C1(m,o,k,l),C2(m,o,k,l),C2s(m,o,k,l),T1max,T2max,"bSSFP signal from "+num2str(ns)+" isochromats, after the "+num2str(l)+" th pulse for fixed $\alpha=$ "+num2str(a(m))+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR(o))+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(k))+" $T_R$","$\frac{t}{T_R}$","");
                ft(m,o,k,l)=ft_;

            end

        end

    end

end