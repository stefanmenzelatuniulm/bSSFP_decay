%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,C1,C2,C2s,T1max,T2max,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+xLabel+" and fitting peaks");

    M=normalize(M,"range");

    coeffs=["T1" "T2" "T2p"];
    fttype = fittype("abs(("+string(C1)+"+"+string(C2)+"*exp(-x/T2p)+"+string(C2s)+"*exp(x/T2p))*exp(-x/T2))",coefficients=coeffs);
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0.001 0.001 0.001],'Upper',[1000 1000 1000],'StartPoint',[13.1 0.6 0.015]);

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
    legend("Normalized signal","Fit with $|C_1(T_1,T_2,\alpha,n)+C_2(T_1,T_2,\alpha,n)e^{-\frac{x}{T_2'}}+C_3(T_1,T_2,\alpha,n)e^{\frac{x}{T_2'}}|e^{-\frac{x}{T_2}}$","interpreter","latex",'fontweight','bold','fontsize',14);

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
    annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', str1, str2, str3],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".png");

    close(fig);

end