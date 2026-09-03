function tau_Q = tauQ_fun(L,PTV,Q,Qin,step)
    
flag_NTCP=0;
    for i=1:size(PTV,1)
        if PTV(i)==1
            if Q(i)<=30/100*Qin(i)
               [rows,cols]=find(L==i);
               for j=1:size(rows)
                   vertex=L(rows(j),:);
                   if Q(vertex)<30/100*Qin(vertex)
                      flag_NTCP=1;
                   end
               end
            end
        end
    end
    
    if flag_NTCP==1
        tau_Q=step;
    else
        tau_Q=0;
    end
end
    


