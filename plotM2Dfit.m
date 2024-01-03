%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,C1,C2,C2s,T1max,T2max,TR,f,ns,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+xLabel+" and fitting");

    
    %abs(x-f*TR)*(1/T2-1/T2p) %x<=TR*f
    %exp(-abs(x-f*TR)*(1/T2+sign(x-TR*f)*1/T2p)) %x>=TR*f

    fitfunction="scale*exp(-abs(x-"+string(f*TR)+")*(1/T2+sign(x-"+string(TR*f)+")*1/T2p))*abs("+string(C1)+"+"+string(C2)+"*exp(-x/T2p)+"+string(C2s)+"*exp(x/T2p))+0*T1";
    coeffs=["T1" "T2p" "T2" "scale"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[1 1 1 0],'Upper',[T1max T2max T2max inf],'StartPoint',[13.1*1000 0.6*1000 0.6*1000 ns]);

    % fitfunction="exp(-x/T2)*abs(C1r+1i*C1i+(C2r+1i*C2i)*exp(-x/T2p)+(C2sr+1i*C2si)*exp(x/T2p))";
    % coeffs=["T2p" "T2" "C1r" "C1i" "C2r" "C2i" "C2sr" "C2si"];
    % options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0 0 0 0 0],'Upper',[T2max T2max inf inf inf inf inf inf],'StartPoint',[0.6*1000 0.6*1000 ns ns ns ns ns ns]);

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
    ylabel("Normalized signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);
    legend("Normalized signal","Fit with $C(N_s)\cdot |C_1(T_1,T_2,\alpha,n)+C_2(T_1,T_2,\alpha,n)e^{-\frac{x}{T_2'}}+C_3(T_1,T_2,\alpha,n)e^{\frac{x}{T_2'}}|e^{-\frac{x}{T_2}}$","interpreter","latex",'fontweight','bold','fontsize',14);

    %legend("Normalized signal","Fit with $|C_1+C_2e^{-\frac{x}{T_2'}}+C_3e^{\frac{x}{T_2'}}|e^{-\frac{x}{T_2}}$","interpreter","latex",'fontweight','bold','fontsize',14);

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{3},coeffvals(3),ci(:,3));
    str4 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(4),ci(:,4));
    % str5 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(5),ci(:,5));
    % str6 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(6),ci(:,6));
    % str7 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(7),ci(:,7));
    % str8 = sprintf('\n %s = %0.3f   (%0.3f   %0.3f)',coeffs{4},coeffvals(8),ci(:,8));
    %annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', str1, str2, str3, str4, str5, str6, str7, str8],'EdgeColor','none',"FitBoxToText","on");
    annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Relaxation constants with 95% confidence bounds: ', str1, str2, str3, str4],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".png");

    close(fig);

end