function [d,mask] = MarrWvlt_Dct(Lx,Ly,M,N,sigma,thr)

% This function is used to construct the dictionary corresponding to the 
% sparse localized component of the decomposition. As shown by the
% instructions document, the dictionary corresponds to the column
% sub-matrix of a circulant matrix. 

% The outputs of the function return the first column of the circulant
% matrix as well as the masking vector, which indicates the columns of the
% circulant matrix that are in the dictionary.

dx = Lx/(N-1);
dy = Ly/(M-1);

flag = 1;
radius = 0;

while flag == 1
    
    r =  (radius*dx)^2 + (radius*dy)^2;
    fval =  (1/(pi*sigma^4)*(1-r/(2*sigma^2)).*exp(-r/(2*sigma^2)));
    
    if abs(fval) <= thr 
        flag = 0;
    else
        radius = radius + 1;
    end 
end

if 2*radius > min(M,N)
    error('The variance of the Marr wavelet is large relative to the sizes of the structure');
end
H = zeros(2*radius+1);

Xc = radius + 1;
Yc = radius + 1;

for i = 1:size(H,1)
    for j = 1:size(H,2)
        r = abs(i - Xc)^2*dx^2 + abs(j - Yc)^2*dy^2; 
        H(i,j) = (1/(pi*sigma^4)*(1-r/(2*sigma^2)).*exp(-r/(2*sigma^2)));
    end
end
% Ind = (abs(H) < thr);
% H(Ind) = 0;
H = H./norm(H,'fro');
d = zeros(M,N);
d( Yc - radius : Yc + radius, Xc - radius : Xc + radius ) = H;
d = d(:)';

mask = zeros(M*N,1);

for i = 1 : N-2*radius

    mask(1+(i-1)*M:i*M-2*radius) = 1;    
    
end

end
