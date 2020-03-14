funcprot(0)

function[Pivots,Pivots_i,a]=U(a_temp,M,N)
   
    k=1
    j_start=1;
    i=1;
    A=a_temp
    pivots_i=1;
    pivots(1,1)=0;
    while(i<=N)
        flag=0
        temp=0;
        for(j=j_start:M)
            if(A(j,i)<>0)
                flag=j;
                temp=i;
                break;
             end
        end
        if(flag<>0 )
            //disp(A,"no swap")
            while(flag<>k)
                    A([flag,flag-1],:)=A([flag-1,flag],:);
                    flag=flag-1;
                    //disp(A,"rowswap")
            end
            k=k+1;
            j_start=k;
            //disp(flag,temp,"flag,temp")
            pivots(pivots_i,1)=flag;
            pivots(pivots_i,2)=temp;
            pivots_i=pivots_i+1;
            
            //making below elements zero
            for(t=flag+1:M)
                //disp("---------->",A(t,temp)/A(flag,temp),"<----")
                 A(t,:)=A(t,:)-(A(t,temp)/A(flag,temp))*A(flag,:);
                 for(zx=1:N)
                     if(A(t,zx)<0.0000000000000001 & A(t,zx)>0)
                         A(t,zx)=0
                     end
                 end 
            end
        end
        //disp(A,"A:")
        i=i+1
     end
     Pivots=pivots
     Pivots_i=pivots_i
     a=A
endfunction

function [Pivots,Pivots_i,a]=R(a_temp,M,N,pivots_i,pivots)
//Conversion of U to R
     A=a_temp;
     if(pivots_i<>0)
         i=pivots_i-1;
         //disp(i,"i:")
         while(i>=1)
             u=pivots(i,1)
             v=pivots(i,2)
             //mprintf("u,v=%d,%d\n",u,v);
             j=u-1;
             while(j>=1)
                 //mprintf("j:%d\n",j); 
                 if(A(j,v)<>0)
                     //disp(A,"before:")
                     //mprintf("->%f\n",(A(j,v)/A(u,v)))
                     A(j,:)=A(j,:)-((A(j,v)/A(u,v))*A(u,:));
                     //disp(A,"after:")
                 end 
                 j=j-1;
             end
             i=i-1;
         end
         for(i=1:pivots_i-1)
             u=pivots(i,1)
             v=pivots(i,2)
             A(u,:)=A(u,:)/A(u,v);
         end
     end
     Pivots=pivots
     Pivots_i=pivots_i
     a=A
endfunction

function space(a_temp,M,N)
    
     [pivots,pivots_i,A]=U(a_temp,M,N)
     disp(A,"U: ")    // U
     mprintf("The pivot columns are : ")
     //disp(pivots)
     for(v=1:(pivots_i-1))
         printf(" Col %d, ",pivots(v,2));
     end
     printf("\n");
     disp("Column space : ")
     mprintf("Basis : {( ")
     for(v=1:(pivots_i-1))
         x=pivots(v,2);
         //disp(a_temp(:,x))
         for(w=1:M)
             mprintf("%f ",a_temp(w,x));
         end
         if(v<>pivots_i-1)
            mprintf(" ),( ");
         end
     end       
     mprintf(")}")
     mprintf("\n");
     mprintf(" Dimension of column space = %d\n ",pivots_i-1);
     mprintf("Rank of the matrix is %d\n",pivots_i-1)
     
     disp("------------------------------------------------------------------")
     disp("ROWSPACE")
     disp("------------------------------------------------------------------")
     
     mprintf("The pivot rows are : ")
     //disp(pivots)
     for(v=1:(pivots_i-1))
         printf(" Row %d, ",pivots(v,1));
     end
     printf("\n");
     
     disp("Row space : ")
     mprintf("Basis : {( ")
     for(v=1:(pivots_i-1))
         x=pivots(v,1);
         //disp(a_temp(:,x))
         for(w=1:N)
             mprintf("%f ",a_temp(x,w));
         end
         if(v<>pivots_i-1)
            mprintf(" ),( ");
         end
     end  
     mprintf(")}")
     mprintf("\n");
     mprintf(" Dimension of row space = %d\n ",pivots_i-1);
     
     disp("----------------------------------------------------------------------\n")
     disp("NULLSPACE")
     disp("------------------------------------------------------------------------")
     
     [pivots,pivots_i,A]=R(A,M,N,pivots_i,pivots)
     disp(A,"Row Reduced form : R =")
     disp("Pivot variables are : ")
     mprintf("Columns : ");
     for(i=1:pivots_i-1)
         mprintf("col %d, ",pivots(i,2));
     end
     mprintf("\n");
     disp("The free variables are : ")
     temp_mat(1)=1;
     for(i=1:N)
         temp_mat(i)=1;
     end
     for(i=1:pivots_i-1)
         temp_mat(pivots(i,2))=0;
     end
     k=1;
     
     free_var_matrix(k)=0;
     for(i=1:N)
         if(temp_mat(i)==1)
             free_var_matrix(k)=i;
             k=k+1;
             mprintf("col %d, ",i);
         end
     end
     mprintf("\n");
     if(k==1)
         mprintf("Basis for NULL SPACE : {( z )}\n ");
         mprintf("Dimension for NULL SPACE=0\n ");  
     
     else    
        mprintf("Basis for NULL SPACE : {( ")
         for(i=1:k-1)
             free_v=free_var_matrix(i);
             //Initialising the basis matrix
             basis_matrix(1)=-1;
             for(j=1:N)
                 basis_matrix(j)=-1;
             end
             for(j=1:k-1)
                 basis_matrix(free_var_matrix(j))=0;
             end
             basis_matrix(free_v)=1;
             //Finding smaller b/w N and M and running loop until then
             N_or_M=N;
             if(M<N)
                 N_or_M=M;
             end
             temp_var=1;
               //loop run till smaller of N and M
             for(j=1:N_or_M)
                basis_elem_val=A(j,free_v);
                while(temp_var<>(N_or_M+1) & basis_matrix(temp_var)<>-1)
                    temp_var=temp_var+1
                end
                // if reached end of basis matrix, break
                if(temp_var==N_or_M+1)
                    break;
                end
                basis_matrix(temp_var)=-1*basis_elem_val;
                temp_var=temp_var+1
              end
              while(temp_var<=N)
                  if(basis_matrix(temp_var)==-1)
                    basis_matrix(temp_var)=0;
                  end
                  temp_var=temp_var+1;
              end
              
              for(j=1:N)
                  mprintf("%f ",basis_matrix(j));
              end
              if(i<>k-1)
                mprintf(" ),( ");
              end
          end  
          
          mprintf(" )}\n");
          mprintf("Dimension of Null space(NULLITY) is : %d\n",k-1);   
      end 
      
     disp("--------------------------------------------------------------------")
     disp("LEFT NULL SPACE")
     disp("--------------------------------------------------------------------")
     
     A=(a_temp)';
     disp("A Transpose = ")
     disp(A)
     Temp=M
     M=N
     N=Temp
     //printf("M:%d N:%d\n",M,N)
     [pivots,pivots_i,A]=U(A,M,N)
     //printf("-<>-")
     [pivots,pivots_i,A]=R(A,M,N,pivots_i,pivots)
     //[pivots,pivots_i,A]=left_null_space(A,M,N,pivots_i,pivots);
     
     disp(A,"R of A transpose")
     disp("Pivot variables are : ")
     mprintf("Columns : ");
     for(i=1:pivots_i-1)
         mprintf("col %d, ",pivots(i,2));
     end
     mprintf("\n");
     disp("The free variables are : ")
        temp_mat(1)=1;
         for(i=1:N)
             temp_mat(i)=1;
         end
         for(i=1:pivots_i-1)
             temp_mat(pivots(i,2))=0;
         end
         k=1;
         
         free_var_matrix(k)=0;
         for(i=1:N)
             if(temp_mat(i)==1)
                 free_var_matrix(k)=i;
                 k=k+1;
                 mprintf("col %d, ",i);
             end
         end
         mprintf("\n");
         if(k==1)
             mprintf("Basis for Left NULL SPACE : {( z )}\n ");
             mprintf("Dimension for Left NULL SPACE=0\n ");  
         
         else    
            mprintf("Basis for LEFT NULL SPACE : {( ")
             for(i=1:k-1)
                 free_v=free_var_matrix(i);
                 //Initialising the basis matrix
                 basis_matrix(1)=-1;
                 for(j=1:N)
                     basis_matrix(j)=-1;
                 end
                 for(j=1:k-1)
                     basis_matrix(free_var_matrix(j))=0;
                 end
                 basis_matrix(free_v)=1;
                 //Finding smaller b/w N and M and running loop until then
                 N_or_M=N;
                 if(M<N)
                     N_or_M=M;
                 end
                 temp_var=1;
                   //loop run till smaller of N and M
                 for(j=1:N_or_M)
                    basis_elem_val=A(j,free_v);
                    while(temp_var<>(N_or_M+1) & basis_matrix(temp_var)<>-1)
                        temp_var=temp_var+1
                    end
                    // if reached end of basis matrix, break
                    if(temp_var==N_or_M+1)
                        break;
                    end
                    basis_matrix(temp_var)=-1*basis_elem_val;
                    temp_var=temp_var+1
                  end
                  while(temp_var<=N)
                      if(basis_matrix(temp_var)==-1)
                        basis_matrix(temp_var)=0;
                      end
                      temp_var=temp_var+1;
                  end
                  
                  for(j=1:N)
                      mprintf("%f ",basis_matrix(j));
                  end
                  if(i<>k-1)
                    mprintf(" ),( ");
                  end
              end  
              
              mprintf(" )}\n");
              mprintf("Dimension of Left Null space is : %d\n",k-1);   
        end 
     
                    
endfunction


function fundamental_spaces()
    disp(" Four Fundamental Subspaces ")
    
    disp("---------------------------------------------------------------------")
    disp("1) COLUMN SPACE")
    disp("--------------------------------------------------------------------")
    disp("Enter the values for the matrix A  ")
    mprintf("Enter values for row size m and column size n \n")
    m=input("m= ")
    n=input("n= ")
    mprintf("Enter values for A :\n")
    for i=1:m
       for j=1:n
           mprintf("Enter value for (%d,%d)",i,j)
           a(i,j)=input("=")
       end
    end
    disp(a,"A:")
    
    space(a,m,n);
   
endfunction

// Main function to be called is fundamental_spaces
fundamental_spaces();


