%Plots dependency of M on 1 variable

function [ftl,ftu]=plotM2Dfit(M,X,plotTitle,xLabel,subfolder)

    disp("Creating 2D plot of M vs "+xLabel+" and fitting peaks");

    M=normalize(M,"range");
    %coeffs=["C1" "C2" "C3" "t2p" "t2"];
    coeffs=["C1r" "C1i" "C2r" "C2i" "C3r" "C3i" "t2p" "t2"];
    %fttype1 = fittype("C1+(C2*exp(-x/t2p)+C3*exp(x/t2p))*exp(-x/t2)",coefficients=coeffs);
    fttype = fittype("abs(C1r+1i*C1i+((C2r+1i*C2i)*exp(-x/t2p)+(C3r+1i*C3i)*exp(x/t2p))*exp(-x/t2))",coefficients=coeffs);
    %options1=fitoptions('Method','NonlinearLeastSquares','Lower',[0 -1000 -1000 0.001 0.001],'Upper',[1 1000 1000 10 10],'StartPoint',[0 1 1 0.6 0.6]);
    options=fitoptions('Method','NonlinearLeastSquares','Lower',[0 -1000000 -1000000 -1000000 -1000000 -1000000 0.001 0.001],'Upper',[1 1000000 1000000 1000000 1000000 1000000 10 10],'StartPoint',[0 1 1 1 1 1 0.6 0.6]);

    [eu,el] = envelope(M,15,"peaks");

    fig=figure('WindowState','maximized');
    
    plot(X,M,"+");
    ax = gca;

    hold on;
    plot(X,eu,"green");

    hold on;
    plot(X,el,"magenta");

    hold on;
    ftu=fit(transpose(X),eu,fttype,options);
    plot(ax,ftu,"r");

    hold on;
    ftl=fit(transpose(X),el,fttype,options);
    plot(ax,ftl,"r");

    ax.FontSize = 14;
    title(plotTitle,"interpreter","latex",'fontweight','bold','fontsize',14);
    xlabel(xLabel,"interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Normalized signal (a. u.)","interpreter","latex",'fontweight','bold','fontsize',14);
    %legend("Signal","Upper envelope of signal","Lower envelope of signal","Upper envelope peaks", "Lower envelope troughs","Fit of upper envelope peaks","Fit of lower envelope troughs");
    legend("Normalized signal","Upper envelope of signal","Lower envelope of signal","Fit of upper envelope","Fit of lower envelope");

    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".fig");
    saveas(fig,pwd+"\Figures\"+subfolder+"\"+strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(strrep(plotTitle," ","_"),".","_"),"$",""),",","_"),"{","_"),"}","_"),"\","_"),"*","_"),"^","_"),"=","_")+".png");

    close(fig);

end