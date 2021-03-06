function z = DTmult(x,dim1,dim2)

% The function computes the two dimensional Discrete Cosine Transform (DCT)
% of the images which are vectorized and saved as columns of x.
% The resulting coefficients are vectorized and stacked in z. 

if size(x,1) ~= (dim1*dim2)
    error('Image dimensions are not correct!');
end

m = size(x,2);

z = zeros(size(x));

for ii = 1 : m
    
    tmp = reshape(x(:,ii),dim1,dim2);
    
    tmp2 = dct2(tmp);
    
    z(:,ii) = tmp2(:); 

end

end