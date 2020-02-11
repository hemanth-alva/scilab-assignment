//Gaussian eliminarion method
format('e',20)
//funcport(0)
A=input("Enter the coeffiecient matrix of (n,n) : ")
b=input("Enter the rigth hand side matrix of (n,1): ")
function [x]=naivegaussianelimination(A,b)
    [n,n1]=size(A);
    [m1,p]=size(b);
    if n~=n1 
        error('Matrix must be a square');
    abort;
    elseif n~=m1
        error('incompactible dimension A & b');
    abort;
end
C=[A b];
//forward elimination
    n=size(A,1);
    for k=1:n-1
        for i=k+1:n
            factor=A(i,k)/A(k,k);
                for j=k+1:n
                    A(i,j)=A(i,j)-factor*A(k,j);
                end
            b(i)=b(i)-factor*b(k);
         end
                 
    end
 //back substitution
       x(n)=b(n)/A(n,n);
       for i=n-1:-1:1
           sum=0;
             for j=i+1:n
                 sum=sum+A(i,j)*x(j); 
             end
   end                  
endfunction
