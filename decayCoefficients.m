%Indices of M: a (1st dimension), TR (2nd dimension), f (3rd dimension), f_eval (4th dimension), pulse (5th dimension)

function [C1 C2 C3 T2 T2p]=decayCoefficients(M)

    dim=size(M);

    for pulse=3:dim(5) %pulse

        for k=1:dim(3) %f
        
            for l=1:dim(1) %a
            
                for m=1:dim(2) %TR


                
                end
            
            end
        
        end
        
    end
    
end