function [DtY] = CTmult(dMtx, maskMtx, Y, M, N)
    
    T = size(Y,2); 

    DtY = [];
    
    for i = 1 : size(maskMtx,2) 
         coli = [dMtx(1,i); flip(dMtx(2:end,i))];
         dum  = real(Circulant_Mult(coli,Y)); 
         Ind  = maskMtx(:,i)==1;
         DtY = [DtY; dum(Ind,:)];
    end
   
end