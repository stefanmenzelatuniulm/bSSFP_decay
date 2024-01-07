clear all;
close all;
clc; 

%-------------SETTINGS-------------

%Number of isochromats
ns=2^15;
            
%Number of -/+ alpha pulses -1 (no sampling after last pulse due to +/-
%alpha/2 tip-back pulse)
n_tot=10;

%Metabolite properties (13C lactate)
w0=300/1000; %in kHz, rotating frame (gamma*B0 is filtered by heterodyne mixing)
FWHM=21/1000; %in kHz %Simulation depends sensitively on FWHM
T1=13.1*1000; %in ms
T2=0.6*1000; %in ms

%Range of flip angles in degree
%Calculation much faster for small cosd(alpha) because steady-state approximation can be applied faster 
amin=55;
amax=55;
na=1;

%Range of repetition times in ms
TRmin=10;
TRmax=10;
nTR=1; 

%Range of times TR*f between initial alpha/2 and first -alpha pulse, if TR
%is the time between -/+ alpha pulses
f=[1/3];

%Calculate signal at time t_eval=f_eval*TR, measured from the end of the
%pulse train. f_eval=0 e.g. calculates the signal directly after the end of
%the pulse train
f_eval_min=0;
f_eval_max=1;
nf_eval=1001; 

%Equilibrium magnetization per isochromat
Meq=1;

%splitfactor loop iterations are used in vectorizedM -> high splitfactor
%causes less RAM usage, but vectorization is not as efficient
splitfactor=128; 

%Recalculate M, or read existing M from save file M.txt?
recalculateM=true;

%Recalculate C1, C2, C2s or read from save file C.txt?
recalculateC=true;

%Steady-state approximation: Set (E1*cosd(alpha))^k = 0 if abs((E1*cos(alpha))^k) < epsilon
%-> Calculations much faster if cosd(alpha) is small
epsilon=0;

%Upper bounds for T1 and T2 (regarding fit, and steady-state approximation)
T1max=20*1000; %in ms
T2max=5*1000; %in ms

%Hyperpolarization factor
hyperpolarization=1;

%-------------END OF SETTINGS-------------

%Flip angles in degree
a=linspace(amin,amax,na);

%Repetition times in ms
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
        M(:,k,:,:,:)=vectorizedM(a,TR(k),w,f,f_eval,n_tot,Meq,T1,T2,hyperpolarization,splitfactor); 
    
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

for k=1:length(f)

    for m=1:na

        for o=1:nTR

            for l=1:n_tot

                M_=permute(M(m,o,k,:,l),[4 1 2 3 5]);

                plotM2Dfit(M_,f_eval,C1(m,o,k,l),C2(m,o,k,l),C2s(m,o,k,l),T1max,T2max,TR,f,ns,"bSSFP signal from "+num2str(ns)+" isochromats, after the "+num2str(l)+" th pulse for fixed $\alpha=$ "+num2str(a(m))+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR(o))+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(k))+" $T_R$","$t$ (ms)","");

            end

        end

    end

end