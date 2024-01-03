%Calculates the coefficients C1, C2, C2s for the decay model
%abs((C1+C2*exp(-t/T2')+C2s*exp(t/T2))*exp(-t/T2)) for all alpha, TR, f and
%pulses. Neglects terms with abs((E1(T1max)*cosd(a))^k)<epsilon -> steady-state
%approximation

function [C1 C2 C2s]=decayCoefficients(as,TR,f,n_tot,epsilon,T1max,Meq)

    disp("Steady-state approximation: Set (E1*cos(alpha))^k = 0 if abs((E1*cos(alpha))^k) < epsilon = "+num2str(epsilon));

    na=length(as);
    nTR=length(TR);
    nf=length(f);

    syms T1 T2 T2p positive; %T1, T2, T2p will later be fit parameters

    n=permute(repmat(linspace(1,n_tot+1,n_tot+1),[na 1 nTR nf]),[1 3 4 2]); %indices: a, TR, f, pulse
    a=permute(repmat(reshape(kron((-1).^(linspace(1,n_tot+1,n_tot+1)+1),as),[na n_tot+1]),[1 1 nTR nf]),[1 3 4 2]); %indices: a, TR, f, pulse
    a(:,:,:,1)=a(:,:,:,1)/2;
    TR=repmat(TR,[na 1 nf n_tot+1]); %indices: a, TR, f, pulse -> TR in ms
    TRf=zeros(na,nTR,nf,n_tot+1);
    TR1mf=zeros(na,nTR,nf,n_tot+1);
    for k=1:nf
        TRf(:,:,k,:)=TR(:,:,k,:)*f(k);
        TR1mf(:,:,k,:)=TR(:,:,k,:)*(1-f(k));
    end
    E1=exp(-TR/T1); %indices: a, TR, f, pulse 
    E2=exp(-TRf*(1/T2-1/T2p)-TR1mf*(1/T2+1/T2p)); %indices: a, TR, f, pulse -> echo bei f*TR: bis dahin T2 decay T2p rephasing, danach T2* decay

    F=Meq*(1-E1); %indices: a, TR, f, pulse
    c1n=1i*sind(a).*E1.*(E1.*cosd(a)).^(n-1)*Meq+1i*sind(a).*E1.*cosd(a).*F.*((1-(E1.*cosd(a)).^(n-1))./(1-E1.*cosd(a)))+1i*sind(a).*F; %indices: a, TR, f, pulse
    c2nm1=E2.*cosd(a/2).^2; %indices: a, TR, f, pulse
    c2snm1=E2.*sind(a/2).^2; %indices: a, TR, f, pulse
    c3nm2mk=0.5*sind(a).^2.*E1.*E2.*(E1.*cosd(a)).^n; %indices: a, TR, f, pulse

    %Steady-state approx. for c1n
    idx=abs((subs(E1,T1,T1max).*cosd(a)).^(n-1))<epsilon;
    c1n(idx)=0;

    %Steady-state approx. for c3nm2mk
    idx=abs((subs(E1,T1,T1max).*cosd(a)).^n)<epsilon;
    c3nm2mk(idx)=0;

    C1=zeros(na,nTR,nf,n_tot+1,"sym"); %indices: a, TR, f, pulse
    C2=zeros(na,nTR,nf,n_tot+1,"sym"); %indices: a, TR, f, pulse
    C2s=zeros(na,nTR,nf,n_tot+1,"sym"); %indices: a, TR, f, pulse

    C1(:,:,:,1)=c1n(:,:,:,1);
    
    for k=2:n_tot+1

        disp("Calculating relaxation coefficients for pulse "+num2str(k-1)+" (counting a/2 preparation pulse)");
    
        C1(:,:,:,k)=c1n(:,:,:,k);
        C2(:,:,:,k)=c2nm1(:,:,:,k).*(C1(:,:,:,k-1)+C2(:,:,:,k-1)+C2s(:,:,:,k-1));
        C2s(:,:,:,k)=c2snm1(:,:,:,k).*(conj(C1(:,:,:,k-1))+conj(C2(:,:,:,k-1))+conj(C2s(:,:,:,k-1)));

        sumC2=zeros(na,nTR,nf);
        
        for l=0:k-3

            if c3nm2mk(:,:,:,l+1)~=0
        
                sumC2=sumC2+c3nm2mk(:,:,:,l+1).*(C1(:,:,:,k-2-l)+C2(:,:,:,k-2-l)+C2s(:,:,:,k-2-l));
            
            end

        end

        C2(:,:,:,k)=C2(:,:,:,k)-sumC2;

        sumC2s=zeros(na,nTR,nf);

        for l=0:k-3

            if c3nm2mk(:,:,:,l+1)~=0
        
                sumC2s=sumC2s+c3nm2mk(:,:,:,l+1).*(conj(C1(:,:,:,k-2-l))+conj(C2(:,:,:,k-2-l))+conj(C2s(:,:,:,k-2-l)));

            end

        end

        C2s(:,:,:,k)=C2s(:,:,:,k)+sumC2s;

    end

    C1=C1(:,:,:,2:end);
    C2=C2(:,:,:,2:end);
    C2s=C2s(:,:,:,2:end);

end