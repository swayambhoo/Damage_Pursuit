function [DtY] = Adjoint_Circulant_Operator(dMtx, maskMtx, Y, M, N)
    
    T = size(Y,2); 

    D1tY = [];
    
    for i = 1 : size(maskMtx,2) 
         coli = [dMtx(1,i); flip(dMtx(2:end,i))];
         dum  = real(multCirculant3(coli,Y)); 
         Ind  = maskMtx(:,i)==1;
         D1tY = [D1tY; dum(Ind,:)];
    end
    
   
   D2tY = [];
   for t = 1 : T
       tmp1 = reshape(Y(:,t), M, N);
       tmp2 = dct2(tmp1);
       D2tY(:,t) = tmp2(:);
   end
   
   DtY = [D1tY; D2tY];
   
end