        function [zHit, zFA, dPrime, c_bias] = GNG_dprime(binsize, go_licks, ngo_licks);
        dPrime = [];
        c_bias = [];
        zHit = [];
        zFA = [];
        
 
            for tt=1:length(go_licks)
                if go_licks(tt)==1;
                    go_licks(tt)=1-1/binsize;
                end
                if ngo_licks(tt)==0;
                    ngo_licks(tt)=1/binsize;
                end
                if go_licks(tt)==0;
                    go_licks(tt)=1/binsize;
                end
                if ngo_licks(tt)==1;
                    ngo_licks(tt)=1-1/binsize;
                end
            end

         
            zHit = norminv(go_licks) ;
            zFA = norminv(ngo_licks) ;
            dPrime = zHit - zFA ; 
           c_bias = (.5*(zFA + zHit)); 



        end
        
        