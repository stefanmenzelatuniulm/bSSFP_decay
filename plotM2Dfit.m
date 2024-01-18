%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,summedTransverseAmplitudes,T1max,T2max,TR,ns,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+strrep(xLabel,"$","")+" and fitting");

    summedTransverseAmplitudes = abs(summedTransverseAmplitudes);

    fitfunction=string(summedTransverseAmplitudes)+"+0*T1+0*T2+0*T2s+0*M_eq";
    coeffs=["T1" "T2" "T2s" "M_eq"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],'Upper',[T1max T2max T2max inf],'StartPoint',[13.1*1000 0.6*1000 1/(pi*21/1000) ns]);

    options2=fitoptions('Method','NonlinearLeastSquares','Lower',[13.1*1000 0.6*1000 1/(pi*21/1000) ns],'Upper',[13.1*1000 0.6*1000 1/(pi*21/1000) ns],'StartPoint',[13.1*1000 0.6*1000 1/(pi*21/1000) ns]);
    
    fttype = fittype(fitfunction,coefficients=coeffs);

    fig=figure('WindowState','maximized');

    X=X*TR; 
    
    plot(X,M,"+");
    ax = gca;

    hold on;
    ft=fit(transpose(X),M,fttype,options);
    ft2=fit(transpose(X),M,fttype,options2);
    pf1 = plot(ax,ft,"r");
    set(pf1,'lineWidth',1);
    hold on;
    pf2 = plot(ax,ft2,"g");
    set(pf2,'lineWidth',1);

    fitfunction3 = "M_eq*exp(-abs(x-5)/T2s)";
    coeffs3=["T2s" "M_eq"];
    options3=fitoptions('Method','NonlinearLeastSquares','Lower',[0 -inf],'Upper',[T2max inf],'StartPoint',[1/(pi*21/1000) ns]);
    fttype3 = fittype(fitfunction3,coefficients=coeffs3);
    ft3=fit(transpose(X),M,fttype3,options3);

    hold on;
    pf3 = plot(ax,ft3,"m");
    set(pf3,'lineWidth',1);

    ax.FontSize = 14;
    title(plotTitle,"interpreter","latex",'fontweight','bold','fontsize',14);
    xlabel(xLabel,"interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Simulated signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);

    latexString = latex(summedTransverseAmplitudes);
    if strlength(latexString)>1100
        latexString = "CharLimit";
    end

    legend("Simulated signal (a. u)","Fit with "+string("$"+latexString+"$"), "'Fit' with known "+"$T_1$"+", $T_2$, $T_2^*$ from simulation settings", "Fit with Mueller model "+string("$"+latex(str2sym("M_eq*exp(-abs(x-5)/T2s)"))+"$"), "interpreter","latex",'fontweight','bold','fontsize',10,"Location","Northwest");
    %legend("Simulated signal (a. u)","Fit", "interpreter","latex",'fontweight','bold','fontsize',14,"Location","Northwest");

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{3},coeffvals(3),ci(:,3));
    str4 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',strrep(coeffs{4},"_",""),coeffvals(4),ci(:,4));
   
    annotation('textbox',[0.521 0.62 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2, str3, str4],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".svg");

    close(fig);

end