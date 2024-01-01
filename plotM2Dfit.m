%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,C1,C2,C2s,T1max,T2max,scale,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+xLabel+" and fitting");

    scalefactor=1000;

    fitfunction=string(scale)+"*abs(("+string(C1)+"+"+string(C2)+"*exp(-x/("+string(scalefactor)+"*T2p))+"+string(C2s)+"*exp(x/("+string(scalefactor)+"*T2p)))*exp(-x/("+string(scalefactor)+"*T2)))";
    fitfunction=strrep(fitfunction,"T1",string(scalefactor)+"*T1");

    if contains(fitfunction,"T1")
        lb=[0 0 0];
        ub=[T1max T2max T2max];
        sp=[13.1 0.6 0.015];
        coeffs=["T1" "T2p" "T2"];
        options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0],'Upper',[T1max/scalefactor T2max/scalefactor T2max/scalefactor],'StartPoint',[13.1/scalefactor 0.6/scalefactor 0.015/scalefactor]);
    else
        lb=[0 0];
        ub=[T2max T2max];
        sp=[0.6 0.015];
        coeffs=["T2p" "T2"];
        options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0],'Upper',[T2max/scalefactor T2max/scalefactor],'StartPoint',[0.6/scalefactor 0.015/scalefactor]);
    end

    fttype = fittype(fitfunction,coefficients=coeffs);

    fig=figure('WindowState','maximized');
    
    plot(X,M,"+");
    ax = gca;

    hold on;
    ft=fit(transpose(X),M,fttype,options);
    plot(ax,ft,"r");

    ax.FontSize = 14;
    title(plotTitle,"interpreter","latex",'fontweight','bold','fontsize',14);
    xlabel(xLabel,"interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Normalized signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);
    legend("Normalized signal","Fit with $N_s|C_1(T_1,T_2,\alpha,n)+C_2(T_1,T_2,\alpha,n)e^{-\frac{x}{T_2'}}+C_3(T_1,T_2,\alpha,n)e^{\frac{x}{T_2'}}|e^{-\frac{x}{T_2}}$","interpreter","latex",'fontweight','bold','fontsize',14);

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1)*scalefactor,ci(:,1)*scalefactor);
    str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2)*scalefactor,ci(:,2)*scalefactor);
    if contains(fitfunction,"T1")
        str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3)*scalefactor,ci(:,3)*scalefactor);
        annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', str1, str2, str3],'EdgeColor','none',"FitBoxToText","on");
    else
        annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', "T1 = NaN", str1, str2],'EdgeColor','none',"FitBoxToText","on");
    end

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".png");

    close(fig);

end