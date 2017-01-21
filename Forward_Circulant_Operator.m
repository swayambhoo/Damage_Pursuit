function [DX] = Forward_Circulant_Operator(dMtx, maskMtx, X, M, N)
    
    T = size(X,2); 
    DX1 = 0;
    count = 0;
    
    for i = 1 : size(maskMtx,2)
    
        Ind  = maskMtx(:,i)==1;
        d = dMtx(:,i);
        r = sum(maskMtx(:,i));
        
        Xtmp = zeros(M*N,T);
        Xtmp(Ind,:) = X( count + 1 : count + r ,:);
        
%         row1 = [d(1) flip(d(end:-1:2))'];
        tmp = real(multCirculant3(d,Xtmp));
        DX1 = DX1 + tmp;
        
        count = count + r;
    
    end
    
   X2 = X(count+1:end,:);
   DX2 = zeros(M*N, T);
   for t = 1 : T
       tmp1 = reshape(X2(:,t), M, N);
       tmp2 = idct2(tmp1);
       DX2(:,t) = tmp2(:);
   end
   
   DX = DX1 + DX2;
   
end