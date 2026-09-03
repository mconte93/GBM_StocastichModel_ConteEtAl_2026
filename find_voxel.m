function [col] = find_voxel(x,y,Pix,h)
%Function that determines in which voxel is located the vertex of
%(x,y)-coordinates.

pos_A = find(Pix(1,:)>=x*2/h-1 & Pix(1,:)<=x*2/h+1);
i=1;
flag=false;

while (i<=size(pos_A,2) & flag==false)
    j=pos_A(i);
    if (Pix(2,j)>=y*2/h-1 & Pix(2,j)<=y*2/h+1)
        col=j;
        flag=true;
    end
    i=i+1;
end
if flag==false 
    col=0;
end
end

