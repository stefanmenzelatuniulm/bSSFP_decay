clear all;
close all;
clc; 

%evtl T1, T2 vernachlässigen gegenüber T2s -> auf inf setzen

%-------------SETTINGS-------------

%Number of isochromats
ns=2^10;
            
%Number of -/+ alpha pulses -1 (no sampling after last pulse due to +/-
%alpha/2 tip-back pulse), not counting a/2 preparation pulse
n_tot=32;

%Assume that steady state is reached after how many pulses
n_steady_state = 12;

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
f=[1/2 1/3];

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
splitfactor=256; 

%Recalculate M, or read existing M from save file M.mat?
recalculateM=true;

%Recalculate transverse Amplitudes or read from save file transverseAmplitudes.mat?
recalculateAmplitudes=true;

%-------------END OF SETTINGS-------------

for g = 1:3

    if g == 1
        oneIso = true;
        finiteIso = false;
        infiniteIso = false;
    elseif g == 2
        oneIso = false;
        finiteIso = true;
        infiniteIso = false;   
        recalculateM = true; 
        recalculateAmplitudes=false;
    else
        oneIso = false;
        finiteIso = false;
        infiniteIso = true;
        recalculateM = true; 
        recalculateAmplitudes=false;
    end
    
    f_eval=linspace(f_eval_min,f_eval_max,nf_eval);
    
    %Lorentzian distribution with FWHM around w0 for Larmor frequencies in Hz
    %(Lorentzian distribution is special case of Student's t distribution with
    %1 degree of freedom) Care: sigma is HWHM
    pd=makedist('tLocationScale','mu',w0,'sigma',FWHM/2,'nu',1);
    w=random(pd,1,ns);
    
    if infiniteIso
        plotHist(w,w0,FWHM,"infiniteIso");
    end
    
    if oneIso
        w = ones(1,ns)*w0;
        plotHist(w,w0,FWHM,"oneIso");
    end
    
    if finiteIso
        w = random(pd,1,128);
        plotHist(w,w0,FWHM,"finiteIso");
    end
    
    %Calculate M Cant fully exploit the vectorization of vectorizedM due to RAM
    %overflow Balance has to be found of CPU time vs RAM usage
    if recalculateM
    
        M=zeros(1,1,length(f),nf_eval,n_tot); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
            
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
    
    for m = 1:length(f)
    
        if recalculateAmplitudes
        
            [transverseAmplitudes, transverseAmplitudesNoIntegration] = sumTransverseAmplitudes(n_tot, a, TR, f(m), n_steady_state, hyperpolarization, g==1, w0);
            save(pwd+"\"+"transverseAmplitudes.mat","transverseAmplitudes","-v7.3");
            save(pwd+"\"+"transverseAmplitudesNoIntegration.mat","transverseAmplitudesNoIntegration","-v7.3");
        
        else
        
            load(pwd+"\"+"transverseAmplitudes.mat");
            load(pwd+"\"+"transverseAmplitudesNoIntegration.mat");
        
        end
        
        disp(" ");
        
        for k=1:n_tot
        
            M_=permute(M(:,:,m,:,k),[4 1 2 3 5]);
    
            if oneIso
                plotM2Dfit(M_,f_eval,transverseAmplitudesNoIntegration(k),T1max,T2max,TR,ns,"bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+" th pulse for fixed $\alpha=$ "+num2str(a)+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$","$t$ (ms)","oneIso\"+num2str(m));
            elseif finiteIso
                plotM2Dfit(M_,f_eval,transverseAmplitudesNoIntegration(k),T1max,T2max,TR,ns,"bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+" th pulse for fixed $\alpha=$ "+num2str(a)+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$","$t$ (ms)","finiteIso\"+num2str(m));
            else           
                plotM2Dfit(M_,f_eval,transverseAmplitudes(k),T1max,T2max,TR,ns,"bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+" th pulse for fixed $\alpha=$ "+num2str(a)+" $^{\circ}$ for fixed $T_R=$ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$","$t$ (ms)","infiniteIso\"+num2str(m));
            end
    
        end
    
    end

end