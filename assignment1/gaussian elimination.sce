//Gaussian eliminarion method
C =zeros(3,4);
//3x4 matrix
for(i=1:3)
    for(j=1:3)
        printf("Enter element a%d%d",i,j);
        C(i,j)=(input(" "));
    end
 end
 for(i=1:3)
     printf("enter b%d",i);
     C(i,4)=(input(" "));
 end
 disp("-------------------------------------------------------------------------------------------------------------------------------------------------------------");
 disp("The augmented matrix is ");
 A=C(1:3,1:3);
 disp(C);
 disp("-------------------------------------------------------------------------------------------------------------------------------------------------------------");

 A=C(1:3,1:3);

 function Gaussian(C)
 
       for(i=1:3)
           disp("-------------------------------------------------------------------------------------------------------------------------------------------------------------");
        if(C(i,i)==0) then
            if( i+1<4 & C(i+1,i)~=0  ) then
                C([i,i+1],:)=C([i+1,i],:);
            
            elseif(i+2<4 & C(i+2,i)~=0 ) then
                 C([i,i+2],:)=C([i+2,i],:);
            end
        end
        if(C(i,i)==0) then
            disp('SINGULAR CASE');break;
        end
        for(j=i+1:3)
            
                C(j,:)=C(j,:)-(C(j,i)/C(i,i))*C(i,:);
        end
        
          end
      
        disp(C);
        
 if(C(3,3)~=0 & C(2,2)~=0 &C(1,1)~=0) then
 x=0;y=0;z=0;
 z=C(3,3)\C(3,4);
 y=(C(2,4) - z*C(2,3))/C(2,2);
 x=(C(1,4) - z*C(1,3) -y*C(1,2))/C(1,1);
 sol=[x,y,z];
 disp(sol,"the solutions are x  y  z repectively ");
end
endfunction
Gaussian(C);
 

