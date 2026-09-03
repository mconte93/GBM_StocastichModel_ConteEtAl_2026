function Q_mesh=matchQ(Q,C,Pix,h) 

%%MATCHQ: Maps fiber data to the computational mesh based on voxel coordinates.

Q_mesh=zeros(3,size(C,2));

for i=1:size(C,2)
    x=C(1,i);
    y=C(2,i);
    [col]=find_voxel(x,y,Pix,h);
    Q_mesh(1:2,i)=C(1:2,i);
    if col==0 
        Q_mesh(3,i)=0;
    else
    Q_mesh(3,i)=Q(3,col);
    end
end

end



