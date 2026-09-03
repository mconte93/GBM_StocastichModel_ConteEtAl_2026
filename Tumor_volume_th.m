function [M,M_th2] = Tumor_volume_th(C,L,PTV,U,th_T2)

%TUMOR_VOLUME_TH: Computes the total tumor mass/volume and the T2-thresholded volume.

s=size(L);     %number of triangles

M_Tumor_T2=0;
M_Tumor=0;

for k=1:s(1)
    K=C(1:2,L(k,:));
    if PTV(L(k,:))==[1 1 1]'
      M_TumorEl=mass_triang(K,U(L(k,:)));
      M_Tumor=M_Tumor+M_TumorEl; 

        if U(L(k,:))>=th_T2
       M_TumorEl_th2=M_TumorEl;
        else
       M_TumorEl_th2=0;
        end
       M_Tumor_T2=M_Tumor_T2+M_TumorEl_th2;  
    end
  
  clear M_TumorEl  
end

M=M_Tumor;
M_th2=M_Tumor_T2;

end