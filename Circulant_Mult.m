function [CX] = Circulant_Mult(c,X)
% This function multiplies circulant matrix whose 
% first column 'c' is with the data matrix 'X'.
        [n,N] = size(X);
        tmp1 = (fft(c));
        tmp2 = bsxfun(@times,tmp1,fft(X));
	CX = ifft(tmp2) ;
end