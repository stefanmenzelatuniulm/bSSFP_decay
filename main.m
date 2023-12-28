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
amin=10;
amax=180;
na=18; %resolution: 10 degree

%Range of repetition times in ms
TRmin=5;
TRmax=20;
nTR=151; %resolution: 1/10 ms

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
recalculateM=true;

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

    M=zeros(na,nTR,length(f),length(f_eval),n_tot); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
    
    for k=1:length(TR)

        M(:,k,:,:,:)=vectorizedM(a,TR(k),w,f,f_eval,n_tot,Meq,T1,T2,splitfactor); 
    
    end

    %Save M and w
    save(pwd+"\"+"M.mat","M","-v7.3");
    save(pwd+"\"+"w.mat","w","-v7.3");

else

    load(pwd+"\"+"M.mat");
    load(pwd+"\"+"w.mat");

end

deleteFigures("Figures");

%Plot histogram of w
plotHist(w,w0,FWHM);

%Plot M(f_eval_fine) for fixed alpha and fixed TR for all pulses separately
for k=1:length(f)

    mkdir("Figures");

    for m=1:length(a)

        for o=1:length(TR)

            for l=1:n_tot

                M_=permute(M(m,o,k,:,l),[4 1 2 3 5]);

                plotM2Dfit(M_,f_eval,"bSSFP signal from "+num2str(ns)+" isochromats, after the "+num2str(l)+" th pulse for fixed $\alpha=$ "+num2str(a(m))+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR(o))+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(k))+" $T_R$","$\frac{t}{T_R}$","");

            end

        end

    end

end