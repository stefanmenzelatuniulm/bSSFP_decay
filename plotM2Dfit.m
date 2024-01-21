%Plots dependency of M on 1 variable

function ft = plotM2Dfit(M, X, summedTransverseAmplitudes, transverseAmplitudesPhaseNoInt, T1, T2, FWHM, TR, ns, plotTitle, xLabel, subfolder, discreteSumInsteadOfIntegral, f)

    disp("Creating 2D plot of M vs "+strrep(xLabel, "$", "")+" and fitting");

    syms w;

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
        
    fitfunction = string(summedTransverseAmplitudes);
    fitfunctionMueller = "C*exp(-abs(x-"+num2str(TR*f)+")/T2s)";
        
    coeffs = ["T1" "T2" "T2s" "M_eq"];
    coeffsMueller = ["T2s" "C"];

    options = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [T1 T2 1/(pi*FWHM) ns], 'Upper', [T1 T2 1/(pi*FWHM) ns], 'StartPoint', [T1 T2 1/(pi*FWHM) ns]);
    optionsT2Dominated = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [realmax T2 1/(pi*FWHM) ns], 'Upper', [realmax T2 1/(pi*FWHM) ns], 'StartPoint', [realmax T2 1/(pi*FWHM) ns]);    
    optionsT2sDominated = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [realmax realmax 1/(pi*FWHM) ns], 'Upper', [realmax realmax 1/(pi*FWHM) ns], 'StartPoint', [realmax realmax 1/(pi*FWHM) ns]);
    optionsMueller = fitoptions('Method', 'NonlinearLeastSquares', 'Lower', [1/(pi*FWHM) 0], 'Upper', [1/(pi*FWHM) inf], 'StartPoint', [1/(pi*FWHM) ns]);
    
    fttype = fittype(fitfunction, coefficients = coeffs);
    fttypeMueller = fittype(fitfunctionMueller, coefficients = coeffsMueller);

    fig = figure('WindowState', 'maximized');

    X = X*TR; 
    plot(X, M, "+");
    ax = gca;

    hold on;
    ft = fit(transpose(X), M, fttype, options);
    pf = plot(ax, ft);
    set(pf, 'lineWidth', 1);
    set(pf, 'color', [0.4660 0.6740 0.1880]);

    hold on;
    ftT2Dominated = fit(transpose(X), M, fttype, optionsT2Dominated);
    pfT2Dominated = plot(ax, ftT2Dominated, "m");
    set(pfT2Dominated, 'lineWidth', 1);

    hold on;
    ftT2sDominated = fit(transpose(X), M, fttype, optionsT2sDominated);
    pfT2sDominated = plot(ax, ftT2sDominated, "c");
    set(pfT2sDominated, 'lineWidth', 1);

    hold on;
    ftMueller = fit(transpose(X), M, fttypeMueller, optionsMueller);
    pfMueller = plot(ax, ftMueller, "r");
    set(pfMueller, 'lineWidth', 1);

    ax.FontSize = 14;
    title(plotTitle, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    xlabel(xLabel, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    ylabel("Simulated signal (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

    shortLatexString = "f \bigl( T_1, \: T_2, \: T_2^*, \: M_{eq} \bigr)";

    legend("Simulated signal", "Exact theoretical signal $$"+shortLatexString+"$$ with known "+"$T_1$"+", $T_2$, $T_2^*$, $M_{eq}$ from simulation settings", "Approximated theoretical signal $$\lim_{T_1 \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2$, $T_2^*$, $M_{eq}$ from simulation settings", "Approximated theoretical signal $$\lim_{T_{1, \: 2} \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2^*$, $M_{eq}$ from simulation settings", "State-of-the-art theoretical signal (Mueller) "+string("$$"+strrep(latex(str2sym(fitfunctionMueller)), "T2s", "T_2^*")+"$$")+" with known $T_2^*$ from simulation settings", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");

    writelines(string(summedTransverseAmplitudes), strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model", WriteMode = "append");
    writelines(string(latex(summedTransverseAmplitudes)), strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model_Latex", WriteMode = "append");
    
    %Difficult to wrap tex/latex strings, therefore plain text
    modelText = text(0.0025, 0, textwrap({strrep(strrep(strrep("f(T1, T2, T2*, M_eq) = "+string(summedTransverseAmplitudes), "M_eq", "M_{eq}"), "T1", "T_1"), "T2", "T_2")}, 2400), 'EdgeColor', 'none', "Color", "black", 'FontSize', 1, 'Units', 'normalized');

    extent = get(modelText).Extent;
    height = extent(4);
    set(modelText, 'Position', get(modelText).Position+[0 height/3 0]);

    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".fig");
    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".svg");

    close(fig);

end