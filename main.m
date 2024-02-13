clear all;
close all;
clc; 

%Einzelner Spin -> kein T2, kein T2* da nichts mit ihm interferieren kann
%Isochromat -> haufen identischer spins -> T2, da die einzelnen Spins durch
%Dipol-Dipol WW interferieren
%Summe aus Isochromaten -> T2 (Interferenz von Spins aus demselben
%Isochromaten) + T2* (Interferenz zwischen verschiedenen Isochromaten, 
%endliche Linienbreite) 
%Tobi fragen!

%-------------SETTINGS-------------

%Number of isochromats
ns = 2^17;
            
%Number of -/+ alpha pulses -1 (no sampling after last pulse due to +/-
%alpha/2 tip-back pulse), not counting a/2 preparation pulse
n_tot = 8;

%Assume that steady state is reached after how many pulses
n_steady_state = 15;

%Metabolite properties (13C lactate)
w0 = 300/1000; %in kHz, rotating frame (gamma*B0 is filtered by heterodyne mixing)
FWHM = 21/1000; %in kHz %Simulation depends sensitively on FWHM
T1 = 13000; %in ms
T2 = 600; %in ms

%Flip angle (deg)
a = 55;

%Repetition time (ms)
TR = 10;

%Range of times TR*f between initial alpha/2 and first -alpha pulse, if TR
%is the time between -/+ alpha pulses
f = [1/3 1/2];

%Calculate signal at time t_eval = f_eval*TR, measured from the end of the
%pulse train. f_eval = 0 e.g. calculates the signal directly after the end of
%the pulse train
f_eval_min = 0;
f_eval_max = 1;
nf_eval = 1001; 

%Equilibrium magnetization per isochromat
Meq = 1;

%Hyperpolarization factor
hyperpolarization = 1;

%splitfactor loop iterations are used in vectorizedM -> high splitfactor
%causes less RAM usage in vectorizedM, but vectorization is not as efficient
splitfactor = 1024; 

%Recalculate M, or read existing M from save file M.mat?
recalculateM = false;

%Recalculate transverse Amplitudes or read from save file transverseAmplitudes.mat?
recalculateAmplitudes = false;

%Sum over known frequency distribution, instead of integrating over the
%whole frequency distribution (only affects fit)
discreteSumInsteadOfIntegral = false;

%Frequency offset due to B0 field inhomogeneity (kHz)
Psi = 0/1000; 

%-------------END OF SETTINGS-------------

f_eval = linspace(f_eval_min, f_eval_max, nf_eval);

%Lorentzian distribution with FWHM around w0 for Larmor frequencies in Hz
%(Lorentzian distribution is special case of Student's t distribution with
%1 degree of freedom) Care: sigma is HWHM
pd = makedist('tLocationScale', 'mu', w0, 'sigma', FWHM/2, 'nu', 1);
w = random(pd, 1, ns);

%Calculate M Cant fully exploit the vectorization of vectorizedM due to RAM
%overflow Balance has to be found of CPU time vs RAM usage
if recalculateM

    M = zeros(1, 1, length(f), nf_eval, n_tot); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
        
    M(:, :, :, :, :) = vectorizedM(a, TR, w, f, f_eval, n_tot, Meq, T1, T2, hyperpolarization, Psi, splitfactor); 
  
    %Save M and w
    save(pwd+"\"+"M.mat", "M", "-v7.3");
    save(pwd+"\"+"w.mat", "w", "-v7.3");

else

    load(pwd+"\"+"M.mat");
    load(pwd+"\"+"w.mat");

end

mkdir("Figures");
deleteFigures("Figures");

%Plot histogram of w
%plotHist(w, w0, FWHM);

for m = 1:length(f)

    if recalculateAmplitudes
    
        [transverseAmplitudes, transverseAmplitudesPhaseNoInt] = sumTransverseAmplitudes(n_tot, a, TR, f(m), n_steady_state, hyperpolarization, true, w0);
        save(pwd+"\"+"transverseAmplitudes.mat", "transverseAmplitudes", "-v7.3");
        save(pwd+"\"+"transverseAmplitudesPhaseNoInt.mat", "transverseAmplitudesPhaseNoInt", "-v7.3");
    
    else
    
        load(pwd+"\"+"transverseAmplitudes.mat");
        load(pwd+"\"+"transverseAmplitudesPhaseNoInt.mat");
    
    end
    
    disp(" ");
    
    for k = 1:n_tot
    
        M_ = permute(M(:, :, m, :, k), [4 1 2 3 5]);

        if k ==  1

            plotM2Dfit(M_, f_eval, transverseAmplitudes(k), transverseAmplitudesPhaseNoInt(k), T1, T2, FWHM, TR, ns, "bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+"st pulse for fixed $\alpha = $ "+num2str(a)+" $^{\circ}$ for fixed $T_R = $ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$", "$t$ (ms)", "", discreteSumInsteadOfIntegral, f(m), Psi);
        
        elseif k ==  2

            plotM2Dfit(M_, f_eval, transverseAmplitudes(k), transverseAmplitudesPhaseNoInt(k), T1, T2, FWHM, TR, ns, "bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+"nd pulse for fixed $\alpha = $ "+num2str(a)+" $^{\circ}$ for fixed $T_R = $ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$", "$t$ (ms)", "", discreteSumInsteadOfIntegral, f(m), Psi);

        elseif k ==  3

            plotM2Dfit(M_, f_eval, transverseAmplitudes(k), transverseAmplitudesPhaseNoInt(k), T1, T2, FWHM, TR, ns, "bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+"rd pulse for fixed $\alpha = $ "+num2str(a)+" $^{\circ}$ for fixed $T_R = $ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$", "$t$ (ms)", "", discreteSumInsteadOfIntegral, f(m), Psi);

        else

            plotM2Dfit(M_, f_eval, transverseAmplitudes(k), transverseAmplitudesPhaseNoInt(k), T1, T2, FWHM, TR, ns, "bSSFP signal from "+num2str(ns)+" isochromats after the "+num2str(k)+"th pulse for fixed $\alpha = $ "+num2str(a)+" $^{\circ}$ for fixed $T_R = $ "+num2str(TR)+" ms for initial $\frac{\alpha}{2}$ pulse spacing "+num2str(f(m))+" $T_R$", "$t$ (ms)", "", discreteSumInsteadOfIntegral, f(m), Psi);
       
        end

    end

end