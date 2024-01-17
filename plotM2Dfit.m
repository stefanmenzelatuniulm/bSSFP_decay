%Plots dependency of M on 1 variable

function ft=plotM2Dfit(M,X,summedTransverseAmplitudes,T1max,T2max,TR,ns,w0,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+strrep(xLabel,"$","")+" and fitting");

    %syms T2s w;

    %gamma_ = FWHM/2;
    %gamma_ = 1/(2*pi*T2s);

    %pw = (1/pi)*gamma_/((w-w0)^2+gamma_^2);

    %summedTransverseAmplitudes = abs(int(summedTransverseAmplitudes*pw, w, -inf, inf));
    summedTransverseAmplitudes = abs(summedTransverseAmplitudes);

    %pd=makedist('tLocationScale','mu',w0,'sigma',FWHM/2,'nu',1);
    %w_values=random(pd,1,ns);

    % load(pwd+"\"+"w.mat");
    % w_values = w;
    % clear w;
    % syms w;

    % summedTransverseAmplitudes_sumOverIsochromats = 0;
    % 
    % for k = 1:ns
    %     summedTransverseAmplitudes_sumOverIsochromats = summedTransverseAmplitudes_sumOverIsochromats + subs(summedTransverseAmplitudes,w,w_values(1,k));
    % end
    % 
    % summedTransverseAmplitudes_sumOverIsochromats = abs(summedTransverseAmplitudes_sumOverIsochromats);

    fitfunction=string(summedTransverseAmplitudes)+"+0*T1+0*T2+0*T2s+0*M_eq";
    coeffs=["T1" "T2" "T2s" "M_eq"];
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],'Upper',[inf inf inf inf],'StartPoint',[13.1*1000 0.6*1000 1/(pi*21/1000) ns]);

    options2=fitoptions('Method','NonlinearLeastSquares','Lower',[13.1*1000 0.6*1000 1/(pi*21/1000) ns],'Upper',[13.1*1000 0.6*1000 1/(pi*21/1000) ns],'StartPoint',[13.1*1000 0.6*1000 1/(pi*21/1000) ns]);
    
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

    legend("Simulated signal (a. u)","Fit", "'Fit' with exact parameters from simulation settings", "interpreter","latex",'fontweight','bold','fontsize',14,"Location","Northwest");
    %legend("Simulated signal (a. u)","Fit", "interpreter","latex",'fontweight','bold','fontsize',14,"Location","Northwest");

    coeffs = coeffnames(ft);
    coeffvals= coeffvalues(ft);
    ci = confint(ft,0.95);
    str1 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{1},coeffvals(1),ci(:,1));
    str2 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{2},coeffvals(2),ci(:,2));
    str3 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',coeffs{3},coeffvals(3),ci(:,3));
    str4 = sprintf('\n %s = %0.9f   (%0.9f   %0.9f)',strrep(coeffs{4},"_",""),coeffvals(4),ci(:,4));
   
    annotation('textbox',[0.521 0.642 0.2 0.2],'String',['Fit coefficients with 95% confidence bounds: ', str1, str2, str3, str4],'EdgeColor','none',"FitBoxToText","on");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".svg");

    close(fig);

end