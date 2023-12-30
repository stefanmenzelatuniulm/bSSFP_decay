%Indices of M: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)

function [C1 C2 C2s]=decayCoefficients(as,TR,f,n_tot,epsilon,T1max,Meq)

    disp("Steady-state approximation: Set (E1*cos(alpha))^k = 0 if abs((E1*cos(alpha))^k) < epsilon = "+num2str(epsilon));

    na=length(as);
    nTR=length(TR);
    nf=length(f);

    n=permute(repmat(linspace(1,n_tot,n_tot),[na 1 nTR nf]),[1 3 4 2]); %indices: a, TR, f, pulse
    a=permute(repmat(reshape(kron((-1).^linspace(1,n_tot,n_tot),as),[na n_tot]),[1 1 nTR nf]),[1 3 4 2]); %indices: a, TR, f, pulse
    TR=repmat(TR,[na 1 nf n_tot]); %indices: a, TR, f, pulse

    syms T1 T2 positive; %T1, T2, T2' will later be fit parameters 

    E1=exp(-TR/T1); %indices: a, TR, f, pulse
    E2=exp(-TR/T2); %indices: a, TR, f, pulse , insgessamt Relaxation mit T2 (nicht T2*)?

    F=Meq*(1-E1); %indices: a, TR, f, pulse

    c1n=1i*sind(a).*E1.*(E1.*cosd(a)).^(n-1)*Meq+1i*sind(a).*E1.*cosd(a).*F.*((1-(E1.*cosd(a)).^(n-1))./(1-E1.*cosd(a)))+1i*sind(a).*F; %indices: a, TR, f, pulse
    c2n=E2.*cosd(a).^2; %indices: a, TR, f, pulse
    c2sn=E2.*sind(a).^2; %indices: a, TR, f, pulse
    c3nm2mk=0.5*sind(a).^2.*E1.*E2.*(E1.*cosd(a)).^n; %indices: a, TR, f, pulse

    %Steady-state approx. for c1n
    idx=abs((subs(E1,T1,T1max).*cosd(a)).^(n-1))<epsilon;
    c1n(idx)=0;

    %Steady-state approx. for c3nm2mk
    idx=abs((subs(E1,T1,T1max).*cosd(a)).^n)<epsilon;
    c3nm2mk(idx)=0;

    disp("Calculating relaxation coefficients for pulse 1");

    C1=zeros(na,nTR,nf,n_tot,"sym"); %indices: a, TR, f, pulse
    C2=zeros(na,nTR,nf,n_tot,"sym"); %indices: a, TR, f, pulse
    C2s=zeros(na,nTR,nf,n_tot,"sym"); %indices: a, TR, f, pulse

    C1(:,:,:,1)=c1n(:,:,:,1);
    
    for k=2:n_tot

        disp("Calculating relaxation coefficients for pulse "+num2str(k));
    
        C1(:,:,:,k)=c1n(:,:,:,k);
        C2(:,:,:,k)=c2n(:,:,:,k).*(C1(:,:,:,k-1)+C2s(:,:,:,k-1)+C2s(:,:,:,k-1));
        C2s(:,:,:,k)=c2sn(:,:,:,k).*(conj(C1(:,:,:,k-1))+conj(C2s(:,:,:,k-1))+conj(C2s(:,:,:,k-1)));

        sumC2=zeros(na,nTR,nf);
        
        parfor l=0:k-3

            if c3nm2mk(:,:,:,l+1)~=0
        
                sumC2=sumC2+c3nm2mk(:,:,:,l+1).*(C1(:,:,:,n_tot-2-(l+1))+C2(:,:,:,n_tot-2-(l+1))+C2s(:,:,:,n_tot-2-(l+1)));
            
            end

        end

        C2(:,:,:,k)=C2(:,:,:,k)-sumC2;

        sumC2s=zeros(na,nTR,nf);

        parfor l=0:k-3

            if c3nm2mk(:,:,:,l+1)~=0
        
                sumC2s=sumC2s+c3nm2mk(:,:,:,l+1).*(conj(C1(:,:,:,n_tot-2-(l+1)))+conj(C2(:,:,:,n_tot-2-(l+1)))+conj(C2s(:,:,:,n_tot-2-(l+1))));

            end

        end

        C2s(:,:,:,k)=C2s(:,:,:,k)+sumC2s;

    end

end