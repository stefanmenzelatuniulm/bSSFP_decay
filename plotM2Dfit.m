%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,summedTransverseAmplitudes,T1max,T2max,TR,ns,plotTitle,xLabel,subfolder)

    T1Scale = 1;
    T2Scale = 1;
    T2pScale = 1;
    M_eqScale = 1;

    disp("Creating 2D plot of M vs "+xLabel+" and fitting");

    summedTransverseAmplitudes = abs(summedTransverseAmplitudes);

    fitfunction=strrep(strrep(strrep(strrep(strrep(string(summedTransverseAmplitudes)+"+0*T1+0*T2+0*T2p+0*M_eq", "T2p", "(T3p/"+num2str(T2pScale)+")"), "T1", "(T1/"+num2str(T1Scale)+")"), "T2", "(T2/"+num2str(T2Scale)+")"), "M_eq", "(M_eq/"+num2str(M_eqScale)+")"), "T3p", "T2p");
    coeffs=["T1" "T2" "T2p" "M_eq"];
    %fitfunction = strrep(fitfunction,"*exp(-abs(x - 10/3)/(T2p/1))","*(-1)*exp(-abs(x - 10/3)/(T2p/1))");
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],'Upper',[T1max*T1Scale T2max*T2Scale T2max*T2pScale 2*ns*M_eqScale],'StartPoint',[13.1*1000*T1Scale 0.6*1000*T2Scale T2pScale/(pi*21/1000) ns*M_eqScale]);

    options2=fitoptions('Method','NonlinearLeastSquares','Lower',[13.1*1000*T1Scale 0.6*1000*T2Scale T2pScale/(pi*21/1000) ns*M_eqScale],'Upper',[13.1*1000*T1Scale 0.6*1000*T2Scale T2pScale/(pi*21/1000) ns*M_eqScale],'StartPoint',[13.1*1000*T1Scale 0.6*1000*T2Scale T2pScale/(pi*21/1000) ns*M_eqScale]);
    
    fttype = fittype(fitfunction,coefficients=coeffs);

    fig=figure('WindowState','maximized');

    X=X*TR; 
    
    plot(X,M,"+");
    ax = gca;

    hold on;
    ft=fit(transpose(X),M,fttype,options);
    ft2=fit(transpose(X),M,fttype,options2);
    plot(ax,ft,"r");
    plot(ax,ft2,"g");

    ax.FontSize = 14;
    title(plotTitle,"interpreter","latex",'fontweight','bold','fontsize',14);
    xlabel(xLabel,"interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Simulated signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);

    legend("Simulated signal (a.u)","Fit","interpreter","latex",'fontweight','bold','fontsize',14);

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1)/T1Scale,ci(:,1)/T1Scale);
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2)/T2Scale,ci(:,2)/T2Scale);
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{3},coeffvals(3)/T2pScale,ci(:,3)/T2pScale);
    str4 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',strrep(coeffs{4},"_",""),coeffvals(4)/M_eqScale,ci(:,4)/M_eqScale);
   
    annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2, str3, str4],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".svg");

    close(fig);

end