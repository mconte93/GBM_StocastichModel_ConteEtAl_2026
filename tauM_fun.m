function tau_M = tauM_fun(PTV,U,step,th)

flag_TCP=0;
    for i=1:size(PTV,1)
        if PTV(i)==1
            if U(i)>th/100
                flag_TCP=1;
            end
        end
    end
    
    if flag_TCP==0
        tau_M=step;
    else
        tau_M=0;
    end
    
end