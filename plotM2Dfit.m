%Plots dependency of M on 1 variable

function plotM2Dfit(M, X, summedTransverseAmplitudes, transverseAmplitudesPhaseNoInt, T1num, T2num, FWHM, TR, ns, plotTitle, xLabel, subfolder, discreteSumInsteadOfIntegral, f, Psinum)

    disp("Creating 2D plot of M vs "+strrep(xLabel, "$", "")+" and fitting");

    syms w T1 T2 T2s Psi M_eq;

    if ~discreteSumInsteadOfIntegral

        summedTransverseAmplitudes = abs(summedTransverseAmplitudes);

    else

        load(pwd+"\"+"w.mat");
        wValues = w;
        clear w;
        syms w;

        summedTransverseAmplitudes_sumOverIsochromats = 0;

        for k = 1:length(wValues)
        
            summedTransverseAmplitudes_sumOverIsochromats = summedTransverseAmplitudes_sumOverIsochromats + subs(transverseAmplitudesPhaseNoInt, w, wValues(k));

        end

        summedTransverseAmplitudes = abs(summedTransverseAmplitudes_sumOverIsochromats);

        %Meq is equilibrium magnetization per isochromats, in this case
        ns = 1;

    end
        
    fitfunctionMueller = "C*(exp(-abs(x-"+num2str(TR*f)+")/T2s)+exp(-abs(x-"+(num2str(TR-TR*f))+")/T2s))";
    coeffsMueller = ["T2s" "C"];
    fttypeMueller = fittype(fitfunctionMueller, coefficients = coeffsMueller);
    optionsMueller = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [1/(pi*FWHM) 0], 'Upper', [1/(pi*FWHM) inf], 'StartPoint', [1/(pi*FWHM) ns]);

    theoreticalSignalExact = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, T1num), T2, T2num), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns);
    theoreticalSignalT2T2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, inf), T2, T2num), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns); 
    theoreticalSignalT1T2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, T1num), T2, inf), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns);
    theoreticalSignalT2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, inf), T2, inf), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns);    
    
    fig = figure('WindowState', 'maximized');

    X = X*TR; 
    plot(X, M, "+");
    ax = gca;

    hold on;
    pExact = fplot(ax, theoreticalSignalExact, [min(X), max(X)]);
    set(pExact, 'lineWidth', 1);
    set(pExact, 'color', [0.4660 0.6740 0.1880]);

    hold on;
    pT2T2s = fplot(ax, theoreticalSignalT2T2s, [min(X), max(X)], "m");
    set(pT2T2s, 'lineWidth', 1);

    hold on;
    pT1T2s = fplot(ax, theoreticalSignalT1T2s, [min(X), max(X)]);
    set(pT1T2s, 'lineWidth', 1);
    set(pT1T2s, 'color', [0.9290 0.6940 0.1250]);

    hold on;
    pT2s = fplot(ax, theoreticalSignalT2s, [min(X), max(X)], "c");
    set(pT2s, 'lineWidth', 1);

    hold on;
    ftMueller = fit(transpose(X), M, fttypeMueller, optionsMueller);
    pfMueller = plot(ax, ftMueller, "r");
    set(pfMueller, 'lineWidth', 1);

    ax.FontSize = 14;
    title(plotTitle, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    xlabel(xLabel, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    ylabel("Simulated signal (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

    shortLatexString = "f \bigl( T_1, \: T_2, \: T_2^*, \: M_{eq}, \: \Psi \bigr)";

    legend("Simulated signal", "Exact theoretical signal $$"+shortLatexString+"$$ with known "+"$T_1$"+", $T_2$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_1 \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_2 \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_1$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_{1, \: 2} \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "State-of-the-art theoretical signal (Mueller) "+string("$$"+strrep(latex(str2sym(fitfunctionMueller)), "T2s", "T_2^*")+"$$")+" with known $T_2^*$ from simulation settings", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");

    writelines(string(summedTransverseAmplitudes), pwd+"\Models\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model", WriteMode = "overwrite");
    writelines(string(latex(summedTransverseAmplitudes)), pwd+"\Models\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model_Latex", WriteMode = "overwrite");
    
    %Difficult to wrap tex/latex strings, therefore plain text
    if strlength(string(summedTransverseAmplitudes))<50000
        modelText = text(0.0025, 0, textwrap({strrep(strrep(strrep("f(T1, T2, T2*, M_eq, \Psi) = "+string(summedTransverseAmplitudes), "M_eq", "M_{eq}"), "T1", "T_1"), "T2", "T_2")}, 2400), 'EdgeColor', 'none', "Color", "black", 'FontSize', 1, 'Units', 'normalized');
    else
        modelText = text(0.0025, 0, "f(T1, T2, T2*, M_eq, \Psi) too long to display", 'EdgeColor', 'none', "Color", "black", 'FontSize', 1, 'Units', 'normalized');
    end

    extent = get(modelText).Extent;
    height = extent(4);
    set(modelText, 'Position', get(modelText).Position+[0 height/3 0]);

    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".fig");
    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".svg");

    close(fig);

end