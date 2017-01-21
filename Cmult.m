function [DX] = Cmult(dMtx, maskMtx, X, M, N)
    
    T = size(X,2); 
    DX = 0;
    count = 0;
    
    for i = 1 : size(maskMtx,2)
    
        Ind  = maskMtx(:,i)==1;
        d = dMtx(:,i);
        r = sum(maskMtx(:,i));
        
        Xtmp = zeros(M*N,T);
        Xtmp(Ind,:) = X( count + 1 : count + r ,:);
        
        tmp = real(Circulant_Mult(d,Xtmp));
        DX = DX + tmp;
        
        count = count + r;
    
    end
   
end