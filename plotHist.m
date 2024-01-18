%Plots histogram for the Larmor frequencies Deltaw of the different
%isochromats in the rotated frame

function plotHist(w,w0,FWHM)  

    [~,ns]=size(w);

    fig=figure('WindowState', 'maximized');
    histogram(w,ns*10);

    xlim([w0-FWHM*20,w0+FWHM*20]);
    title("Histogram of $\Delta\omega$ of "+num2str(ns)+" different isochromats","interpreter","latex",'fontweight','bold','fontsize',16);
    xlabel("$\Delta\omega$ (kHz)","interpreter","latex",'fontweight','bold','fontsize',14);
    ylabel("Occurrences","interpreter","latex",'fontweight','bold','fontsize',14);
    
    ax = gca;
    ax.FontSize = 14; 

    saveas(fig,pwd+"/Figures/"+"Histw"+num2str(ns)+"isochromats.fig");
    saveas(fig,pwd+"/Figures/"+"Histw"+num2str(ns)+"isochromats.svg");
    
    close(fig);

end
