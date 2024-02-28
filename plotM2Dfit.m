%Plots dependency of M on 1 variable

function plotM2Dfit(M, X, summedTransverseAmplitudes, transverseAmplitudesPhaseNoInt, T1num, T2num, FWHM, TR, ns, plotTitle, xLabel, subfolder, discreteSumInsteadOfIntegral, f, Psinum)

    disp("Creating 2D plot of M vs "+strrep(xLabel, "$", "")+" and fitting");

    syms w T1 T2 Psi M_eq;
    T2s = sym("T2s");

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
    theoreticalSignalT2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, inf), T2, inf), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns);  
    shortLatexString = "f \bigl( T_1, \: T_2, \: T_2^*, \: M_{eq}, \: \Psi \bigr)";

    % if ~fitOnly

        % theoreticalSignalT2T2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, inf), T2, T2num), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns); 
        % theoreticalSignalT1T2s = subs(subs(subs(subs(subs(summedTransverseAmplitudes, T1, T1num), T2, inf), T2s, 1/(pi*FWHM)), Psi, Psinum), M_eq, ns);  
        % 
        % fig = figure('WindowState', 'maximized');
        % 
        % X = X*TR; 
        % plot(X, M, "+");
        % ax = gca;
        % 
        % hold on;
        % pExact = fplot(ax, theoreticalSignalExact, [min(X), max(X)]);
        % set(pExact, 'lineWidth', 1);
        % set(pExact, 'color', [0.4660 0.6740 0.1880]);
        % 
        % hold on;
        % pT2T2s = fplot(ax, theoreticalSignalT2T2s, [min(X), max(X)], "m");
        % set(pT2T2s, 'lineWidth', 1);
        % 
        % hold on;
        % pT1T2s = fplot(ax, theoreticalSignalT1T2s, [min(X), max(X)]);
        % set(pT1T2s, 'lineWidth', 1);
        % set(pT1T2s, 'color', [0.9290 0.6940 0.1250]);
        % 
        % hold on;
        % pT2s = fplot(ax, theoreticalSignalT2s, [min(X), max(X)], "c");
        % set(pT2s, 'lineWidth', 1);
        % 
        % hold on;
        % ftMueller = fit(transpose(X), M, fttypeMueller, optionsMueller);
        % pfMueller = plot(ax, ftMueller, "r");
        % set(pfMueller, 'lineWidth', 1);
        % 
        % ax.FontSize = 14;
        % title(plotTitle, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
        % xlabel(xLabel, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
        % ylabel("Simulated signal (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
        % 
        % legend("Simulated signal", "Exact theoretical signal $$"+shortLatexString+"$$ with known "+"$T_1$"+", $T_2$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_1 \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_2 \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_1$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_{1, \: 2} \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "State-of-the-art theoretical signal (Mueller) "+string("$$"+strrep(latex(str2sym(fitfunctionMueller)), "T2s", "T_2^*")+"$$")+" with known $T_2^*$ from simulation settings", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");
        % 
    % else

    fitfunction=string(summedTransverseAmplitudes)+"+0*T1+0*T2+0*T2s+0*M_eq+0*Psi";
    coeffs=["T1" "T2" "T2s" "M_eq" "Psi"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0 0],'Upper',[inf inf inf inf inf],'StartPoint',[13.1*1000 0.6*1000 1/(pi*21/1000) ns 0]);
    fttype = fittype(fitfunction,coefficients=coeffs);

    fitfunction2 = string(subs(subs(summedTransverseAmplitudes, T1, inf), T2, inf))+"+0*T2s+0*M_eq+0*Psi"; 
    coeffs2=["T2s" "M_eq" "Psi"];
    options2=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[inf inf inf],'StartPoint',[1/(pi*21/1000) ns 0]);
    fttype2 = fittype(fitfunction2,coefficients=coeffs2);
    
    fig = figure('WindowState', 'maximized');

    X = X*TR; 
    plot(X, M, "+");
    ax = gca;

    hold on;
    pExact = fplot(ax, theoreticalSignalExact, [min(X), max(X)]);
    set(pExact, 'lineWidth', 1);
    set(pExact, 'color', "r");

    hold on;
    pT2s = fplot(ax, theoreticalSignalT2s, [min(X), max(X)], "c");
    set(pT2s, 'lineWidth', 1);

    hold on;
    ft=fit(transpose(X),M,fttype,options);
    pf1 = plot(ax,ft,"r");
    set(pf1, 'lineWidth', 1);
    set(pf1, 'color', [0.4660 0.6740 0.1880]); 

    hold on;
    ft2=fit(transpose(X),M,fttype2,options2);
    pf2 = plot(ax,ft2,"r");
    set(pf2, 'lineWidth', 1);
    set(pf2, 'color', [0.9290 0.6940 0.1250]); 

    hold on;
    ftMueller = fit(transpose(X), M, fttypeMueller, optionsMueller);
    pfMueller = plot(ax, ftMueller, "m");
    set(pfMueller, 'lineWidth', 1);

    ax.FontSize = 14;
    title(plotTitle, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    xlabel(xLabel, "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);
    ylabel("Signal (a. u.)", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 14);

    legend("Simulated signal", "Exact theoretical signal $$"+shortLatexString+"$$ with known "+"$T_1$"+", $T_2$, $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Approximated theoretical signal $$\lim_{T_{1, \: 2} \: \to \: \infty} \: "+shortLatexString+"$$ with known $T_2^*$, $M_{eq}$, $\Psi$ from simulation settings", "Fit with $$"+shortLatexString+"$$ with unknown "+"$T_1$"+", $T_2$, $T_2^*$, $M_{eq}$, $\Psi$", "Fit with $$\lim_{T_{1, \: 2} \: \to \: \infty} \: "+shortLatexString+"$$ with unknown $T_2^*$, $M_{eq}$, $\Psi$", "Fit with state-of-the-art theoretical signal $$"+"f_{Mueller}(T_2^*,C)="+latex(str2sym(fitfunctionMueller))+"$$ with unknown $T_2^*$, $C$", "interpreter", "latex", 'fontweight', 'bold', 'fontsize', 10, "Location", "Northwest");

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{3},coeffvals(3),ci(:,3));
    str4 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',strrep(coeffs{4},"_",""),coeffvals(4),ci(:,4));
    str5 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{5},coeffvals(5),ci(:,5));

    annotation('textbox',[0.53 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2, str3, str4, str5],'EdgeColor','none',"FitBoxToText","on", "Color",[0.4660 0.6740 0.1880],"FontSize",4);

    coeffs = coeffnames(ft2);
    coeffvals= coeffvalues(ft2);
    ci = confint(ft2,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',strrep(coeffs{2},"_",""),coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{3},coeffvals(3),ci(:,3));

    annotation('textbox',[0.65 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2, str3],'EdgeColor','none',"FitBoxToText","on", "Color",[0.9290 0.6940 0.1250],"FontSize",4);

    coeffs = coeffnames(ftMueller);
    coeffvals= coeffvalues(ftMueller);
    ci = confint(ftMueller,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2));

    annotation('textbox',[0.77 0.69 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2],'EdgeColor','none',"FitBoxToText","on", "Color","m","FontSize",4);

    writelines(string(summedTransverseAmplitudes), pwd+"\Models\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model", WriteMode = "overwrite");
    writelines(string(latex(summedTransverseAmplitudes)), pwd+"\Models\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+"_Model_Latex", WriteMode = "overwrite");
    
    % %Difficult to wrap tex/latex strings, therefore plain text
    % if strlength(string(summedTransverseAmplitudes))<50000
    %     modelText = text(0.0025, 0, textwrap({strrep(strrep(strrep("f(T1, T2, T2*, M_eq, \Psi) = "+string(summedTransverseAmplitudes), "M_eq", "M_{eq}"), "T1", "T_1"), "T2", "T_2")}, 2400), 'EdgeColor', 'none', "Color", "black", 'FontSize', 1, 'Units', 'normalized');
    % else
    %     modelText = text(0.0025, 0, "f(T1, T2, T2*, M_eq, \Psi) too long to display", 'EdgeColor', 'none', "Color", "black", 'FontSize', 1, 'Units', 'normalized');
    % end
    % 
    % extent = get(modelText).Extent;
    % height = extent(4);
    % set(modelText, 'Position', get(modelText).Position+[0 height/3 0]);

    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".fig");
    saveas(fig, pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle, " ", "_"), ".", "_"), "$", ""), ", ", "_"), "{", "_"), "}", "_"), "\", "_"), "*", "_"), "^", "_"), " = ", "_")+".svg");

    close(fig);

end