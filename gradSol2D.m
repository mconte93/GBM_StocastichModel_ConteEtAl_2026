function GradU = gradSol2D(C,L,U)

%%%%%%%%%%%%%%%%%%%%%%%
% GRADU: Calculation of the gradient of the solution as an average of the different triangle on each node
%%%%%%%%%%%%%%%%%%%%%%

GradU=zeros(size(C));
    for i=1:size(C,2)
        [row,col]=find(L==i);   %Finding the elements to which the node i belongs
        gU=[0;0];

        for k=1:length(row)
            K=C(:,L(row(k),:));
            gradU=DivU(K,U(L(row(k),:)));   %Approximation of the gradient on each element
            gU=gU+gradU(:,col(k)); 
        end
        
        GradU(:,i)=gU./length(row);    %mean of the founded values
        clear row col gU
    end

end