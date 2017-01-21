function y = soft(x,T)

% Soft-thresholding function: y = sign(x) * max((abs(x) - T),0).

if sum(abs(T(:)))==0
   y = x;
else
   y1 = max(abs(x) - T, 0);
   y = ( y1./(y1+T) ).* x;
end


