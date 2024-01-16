%falls alle pulse +a -> kein min bei t=1/2 bei z.B. n=3

clear all;
close all;
clc; 

%-------------SETTINGS-------------

%Number of isochromats
ns=2^15;
            
%Number of -/+ alpha pulses -1 (no sampling after last pulse due to +/-
%alpha/2 tip-back pulse), not counting a/2 preparation pulse
n_tot=3;

%Assume that steady state is reached after how many pulses
n_steady_state = 20;

%Metabolite properties (13C lactate)
w0=300/1000; %in kHz, rotating frame (gamma*B0 is filtered by heterodyne mixing)
FWHM=21/1000; %in kHz %Simulation depends sensitively on FWHM
T1=13.1*1000; %in ms
T2=0.6*1000; %in ms

%Flip angle (deg)
a = 55;

%Repetition time (ms)
TR = 10;

%Range of times TR*f between initial alpha/2 and first -alpha pulse, if TR
%is the time between -/+ alpha pulses
f=1/2;

%Calculate signal at time t_eval=f_eval*TR, measured from the end of the
%pulse train. f_eval=0 e.g. calculates the signal directly after the end of
%the pulse train
f_eval_min=0;
f_eval_max=1;
nf_eval=1001; 

%Equilibrium magnetization per isochromat
Meq=1;

%Hyperpolarization factor
hyperpolarization=1;

%Upper fit bounds for T1 and T2 
T1max=inf; %in ms
T2max=inf; %in ms

%splitfactor loop iterations are used in vectorizedM -> high splitfactor
%causes less RAM usage in vectorizedM, but vectorization is not as efficient
splitfactor=128; 

%Recalculate M, or read existing M from save file M.mat?
recalculateM=true;

%Recalculate transverse Amplitudes or read from save file transverseAmplitudes.mat?
recalculateAmplitudes=true;

%-------------END OF SETTINGS-------------

f_eval=linspace(f_eval_min,f_eval_max,nf_eval);

%Lorentzian distribution with FWHM around w0 for Larmor frequencies in Hz
%(Lorentzian distribution is special case of Student's t distribution with
%1 degree of freedom) Care: sigma is HWHM
pd=makedist('tLocationScale','mu',w0,'sigma',FWHM/2,'nu',1);
w=random(pd,1,ns);

%Calculate M Cant fully exploit the vectorization of vectorizedM due to RAM
%overflow Balance has to be found of CPU time vs RAM usage
if recalculateM

    M=zeros(1,1,1,nf_eval,n_tot); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
        
    M(:,:,:,:,:)=vectorizedM(a,TR,w,f,f_eval,n_tot,Meq,T1,T2,hyperpolarization,splitfactor); 
  
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
%plotHist(w,w0,FWHM);

if recalculateAmplitudes

    transverseAmplitudes = sumTransverseAmplitudes(n_tot, a, TR, f, n_steady_state, hyperpolarization, true);
    save(pwd+"\"+"transverseAmplitudes.mat","transverseAmplitudes","-v7.3");

else

    load(pwd+"\"+"transverseAmplitudes.mat");

end

disp(" ");

for k=1:n_tot

    M_=permute(M(:,:,:,:,k),[4 1 2 3 5]);

    plotM2Dfit(M_,f_eval,transverseAmplitudes(k),T1max,T2max,TR,ns,"bSSFP signal from "+num2str(ns)+" isochromats, after the "+num2str(k)+" th pulse for fixed $\alpha=$ "+num2str(a)+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f)+" $T_R$","$t$ (ms)","");

end