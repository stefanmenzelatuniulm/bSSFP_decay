%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,summedTransverseAmplitudes,T1max,T2max,TR,ns,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+xLabel+" and fitting");

    syms T1 T2p T2 M_eq x;
    summedTransverseAmplitudes = summedTransverseAmplitudes*exp(-x/T2)*exp(-x/T2p);

    fitfunction=string(summedTransverseAmplitudes)+"+0*T1+0*T2+0*T2p+0*M_eq";
    coeffs=["T1" "T2p" "T2" "M_eq"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[1 1 1 0],'Upper',[T1max T2max T2max inf],'StartPoint',[13.1*1000 0.6*1000 0.6*1000 ns]);

    fttype = fittype(fitfunction,coefficients=coeffs);

    fig=figure('WindowState','maximized');

    X=X*TR; 
    
    plot(X,M,"+");
    ax = gca;

    hold on;
    ft=fit(transpose(X),M,fttype,options);
    plot(ax,ft,"r");

    ax.FontSize = 14;
    title(plotTitle,"interpreter","latex",'fontweight','bold','fontsize',14);
    xlabel(xLabel,"interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Simulated signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);

    legend("Simulated signal (a.u)","Fit","interpreter","latex",'fontweight','bold','fontsize',14);

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
    str4 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(4),ci(:,4));

    annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', str1, str2, str3, str4],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".svg");

    close(fig);

end