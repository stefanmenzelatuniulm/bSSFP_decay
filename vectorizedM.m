%Calculates magnetization for flip angles a, repetition times TR, Larmor
%frequencies w, initial a/2 pulse followed after t=f*TR by n_tot -/+ alpha
%pulses. Evaluates the signal at t=f_eval*TR, with t measured after the end
%of the pulse train. f*TR is the spacing between the initial a/2 pulse and
%the first -a pulse. Returns M_tot=M(a,TR,f,f_eval,pulse). Will calculate
%M_tot for every isochromat and then sum the result over the isochromats

function M_tot = vectorizedM(a,TR,w,f,f_eval,n_tot,Meq,T1,T2,splitfactor)

    %Convert TR to s
    TR=TR/1000;

    %Number of flip angles, repetition times, isochromats, alpha/2 pulse
    %spacings, evaluation times
    [~,na]=size(a);
    [~,nTR]=size(TR);
    [~,ns]=size(w);
    [~,nf]=size(f);
    [~,nf_eval]=size(f_eval);

    %How many isochromats are considered per loop iteration (only affects
    %balance between RAM usage and computation speed)
    splitfactor=max(splitfactor,1);
    Nwsplit=max(round(ns/splitfactor),1);
    ii=(1:numel(w))';
    wsplit=transpose(accumarray([rem(ii-1,Nwsplit)+1,ceil(ii/Nwsplit)],w,[],[],nan)); %w split into Nwsplit columns, possibly with NaN at the end 
    [wsplitsize,ns]=size(wsplit);

    %Unity vector in z direction
    ez=[0; 0; 1];

    %f*TR after an initial alpha/2 pulse, a pulse train with n_tot
    %alternating -/+ alpha pulses with repetition time TR occurs.
    %M_n(:,1,k,l,m,n,o,p) will contain the transverse magnetizations of
    %k-th flip angle, the l-th repetition time, the m-th a/2 pulse spacing,
    %the n-th evaluation time (0<=fraction<=1 of TR) and the o-th pulse
    M_tot=zeros(3,1,na,nTR,nf,nf_eval,n_tot);

    %alpha (in degree) rotation matrix around x axis
    Rx=@(alpha) [1 0 0; 0 cosd(alpha) sind(alpha); 0 -sind(alpha) cosd(alpha)];

    %alpha rotation matrix
    Rx_a=arrayfun(Rx,a,'UniformOutput',false);
    Rx_a=cell2mat(reshape(Rx_a,[1 1 na])); %indices R 3x3, a
    Rx_a=repmat(Rx_a,[1 1 1 nTR Nwsplit nf nf_eval]); %indices: relaxation 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    %-alpha rotation  matrix
    Rx_ma=arrayfun(Rx,-a,'UniformOutput',false);
    Rx_ma=cell2mat(reshape(Rx_ma,[1 1 na])); %indices R 3x3, a
    Rx_ma=repmat(Rx_ma,[1 1 1 nTR Nwsplit nf nf_eval]); %indices: relaxation 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    %alpha/2 rotation matrix
    Rx_a2=arrayfun(Rx,a/2,'UniformOutput',false);
    Rx_a2=cell2mat(reshape(Rx_a2,[1 1 na])); %indices R 3x3, a
    Rx_a2=repmat(Rx_a2,[1 1 1 nTR Nwsplit nf nf_eval]); %indices: relaxation 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    %Relaxation/Dephasing matrices T1 relaxation terms
    E1r=@(T) [0 0 0; 0 0 0; 0 0 1-exp(-T/T1)];

    E1r_TR=cell2mat(reshape(arrayfun(E1r,TR,'UniformOutput',false),[1 1 nTR]));
    E1r_TR=repmat(E1r_TR,[1 1 1 na Nwsplit nf nf_eval]); %indices: E1r_TR 3x3, TR (3rd dimension), a (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)
    E1r_TR=permute(E1r_TR,[1 2 4 3 5 6 7]); %indices: E1r_TR 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)
    
    TRf=reshape(kron(f,TR),[nTR nf]);
    E1r_TRf=cell2mat(reshape(arrayfun(E1r,TRf,'UniformOutput',false),[1 1 nTR nf]));
    E1r_TRf=repmat(E1r_TRf,[1 1 1 1 na Nwsplit nf_eval]); %indices: E1r_TRf 3x3, TR (3rd dimension), f (4th dimension), a (5th dimension), w (6th dimension), f_eval (7th dimension)
    E1r_TRf=permute(E1r_TRf,[1 2 5 3 6 4 7]); %indices: E1r_TRf 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    TRf_eval=reshape(kron(f_eval,TR),[nTR nf_eval]);
    E1r_TRf_eval=cell2mat(reshape(arrayfun(E1r,TRf_eval,'UniformOutput',false),[1 1 nTR nf_eval]));
    E1r_TRf_eval=repmat(E1r_TRf_eval,[1 1 1 1 na Nwsplit nf]); %indices: E1r_TRf_eval 3x3, TR (3rd dimension), f_eval (4th dimension), a (5th dimension), w (6th dimension), f (7th dimension)
    E1r_TRf_eval=permute(E1r_TRf_eval,[1 2 5 3 6 7 4]); %indices: E1r_TRf_eval 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    relaxation=@(T) [exp(-T/T2) exp(-T/T2) 0; exp(-T/T2) exp(-T/T2) 0; 0 0 exp(-T/T1)];
    
    relaxation_TR=cell2mat(reshape(arrayfun(relaxation,TR,'UniformOutput',false),[1 1 nTR])); %indices: relaxation_TR 3x3, TR (3rd dimension)
    relaxation_TR=repmat(relaxation_TR,[1 1 1 Nwsplit na nf nf_eval]); %indices: relaxation_TR 3x3, TR (3rd dimension), w (4th dimension), a (5th dimension), f (6th dimension), f_eval (7th dimension)
    relaxation_TR=permute(relaxation_TR,[1 2 5 3 4 6 7]); %indices: relaxation_TR 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    relaxation_TRf=cell2mat(reshape(arrayfun(relaxation,TRf,'UniformOutput',false),[1 1 nTR nf]));
    relaxation_TRf=repmat(relaxation_TRf,[1 1 1 1 na Nwsplit nf_eval]); %indices: relaxation_TRf 3x3, TR (3rd dimension), f (4th dimension), a (5th dimension), w (6th dimension), f_eval (7th dimension)
    relaxation_TRf=permute(relaxation_TRf,[1 2 5 3 6 4 7]); %indices:relaxation_TRf 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)

    relaxation_TRf_eval=cell2mat(reshape(arrayfun(relaxation,TRf_eval,'UniformOutput',false),[1 1 nTR nf_eval]));
    relaxation_TRf_eval=repmat(relaxation_TRf_eval,[1 1 1 1 na Nwsplit nf n_tot]); %indices: relaxation_TRf_eval 3x3, TR (3rd dimension), f_eval (4th dimension), a (5th dimension), w (6th dimension), f (7th dimension), pulse (8th dimension)
    relaxation_TRf_eval=permute(relaxation_TRf_eval,[1 2 5 3 6 7 4 8]); %indices: relaxation_TRf_eval 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension), pulse (8th dimension)

    %parfor possible here because the order of summation does not matter
    parfor k=1:wsplitsize

        w=wsplit(k,:);        

        disp("Calculating M for "+num2str(ns)+" isochromats, starting at the "+num2str((k-1)*Nwsplit+1)+" th isochromat");
        disp("Loop iteration: "+num2str(k)+"/"+num2str(wsplitsize)+", f_eval="+num2str(f_eval));
    
        disp("  Relaxation + dephasing matrices for "+num2str(ns)+" isochromats, starting at the "+num2str((k-1)*Nwsplit+1)+" th isochromat");
        dephasing=@(b) [cosd(b) sind(b) 0; -sind(b) cosd(b) 0; 0 0 1];
        wT_=kron(360*w,TR); %Care: dont apply mod(wT_,360) because wTf will use wT_*f
        wT=mod(reshape(wT_,[nTR ns]),360);
        wTf=mod(reshape(kron(f,wT_),[nTR ns nf]),360);
        wTf_eval=mod(reshape(kron(f_eval,wT_),[nTR ns nf_eval]),360);
    
        dephasing_TR=repmat(cell2mat(reshape(arrayfun(dephasing,wT,'UniformOutput',false),[1 1 nTR ns])),[1 1 1 1 na nf nf_eval]);  %indices: dephasing 3x3, TR (3rd dimension), w (4th dimension), a (5th dimension), f (6th dimension), f_eval (7th dimension)
        dephasing_TR=permute(dephasing_TR,[1 2 5 3 4 6 7]);
        D_TR=relaxation_TR.*dephasing_TR; %indices: relaxation*dephasing 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)
    
        dephasing_f=repmat(cell2mat(reshape(arrayfun(dephasing,wTf,'UniformOutput',false),[1 1 nTR ns nf])),[1 1 1 1 1 na nf_eval]);  %indices: dephasing 3x3, TR (3rd dimension), w (4th dimension), f (5th dimension), a (6th dimension), f_eval (7th dimension)
        dephasing_f=permute(dephasing_f,[1 2 6 3 4 5 7]); 
        D_TRf=relaxation_TRf.*dephasing_f; %indices: relaxation*dephasing 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)
        
        dephasing_f_eval=repmat(cell2mat(reshape(arrayfun(dephasing,wTf_eval,'UniformOutput',false),[1 1 nTR ns nf_eval])),[1 1 1 1 1 na nf n_tot]);  %indices: dephasing 3x3, TR (3rd dimension), w (4th dimension), f_eval (5th dimension), a (6th dimension), f (7th dimension), pulse (8th dimension)
        dephasing_f_eval=permute(dephasing_f_eval,[1 2 6 3 4 7 5 8]); 
        D_TRf_eval=relaxation_TRf_eval.*dephasing_f_eval; %indices: relaxation*dephasing 3x3, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension), pulse (8th dimension)
               
        %Magnetization after the initial alpha/2 preparation pulse
        disp("  M0, M1 for "+num2str(ns)+" isochromats, starting at the "+num2str((k-1)*Nwsplit+1)+" th isochromat");
        Ma_0=pagemtimes(Rx_a2,ez)*Meq; %indices: Ma_0 3x1, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension)
        
        %Magnetization after the first -alpha pulse
        Ma_1=pagemtimes(Rx_ma,pagemtimes(D_TRf,Ma_0)+Meq*pagemtimes(E1r_TRf,ez));
 
        %Initialize magnetization directly after the nth pulse
        Ma_n=zeros(3,1,na,nTR,ns,nf,nf_eval,n_tot); %indices: Ma_n 3x1, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension), pulse (8th dimension)
        Ma_n(:,:,:,:,:,:,:,1)=Ma_1; %indices: Ma_n 3x1, a (3rd dimension), TR (4th dimension), w (5th dimension), f (6th dimension), f_eval (7th dimension), pulse (8th dimension)
 
        %Calculate Ma_n with n>=2
        for n=2:n_tot
            
            if mod(int8(n),2)==0

                Ma_n(:,:,:,:,:,:,:,n)=pagemtimes(Rx_a,pagemtimes(D_TR,Ma_n(:,:,:,:,:,:,:,n-1))+Meq*pagemtimes(E1r_TR,ez));

            else

                Ma_n(:,:,:,:,:,:,:,n)=pagemtimes(Rx_ma,pagemtimes(D_TR,Ma_n(:,:,:,:,:,:,:,n-1))+Meq*pagemtimes(E1r_TR,ez));

            end

            % Tip-back pulse after the end of the pulse train is not
            % calculated. Must be Rx(-alpha/2) if n_tot is even, and
            % Rx(alpha/2) if n_tot is odd

        end

        %Sum over isochromats and evaluate at t=f_eval*TR
        disp("  Evaluating at t=f_eval*TR and summing the results for "+num2str(ns)+" isochromats, starting at the "+num2str((k-1)*Nwsplit+1)+" th isochromat");
        M_tot=M_tot+permute(sum(pagemtimes(D_TRf_eval,Ma_n)+Meq*pagemtimes(E1r_TRf_eval,ez),5,"omitnan"),[1 2 3 4 6 7 8 5]); 

    end

    M_tot=sqrt(M_tot(1,1,:,:,:,:,:).^2+M_tot(2,1,:,:,:,:,:).^2);
    M_tot=permute(M_tot,[3 4 5 6 7 1 2]); %indices: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)
        
end