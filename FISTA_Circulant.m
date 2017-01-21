function [Beta1,Beta2,f,time] = FISTA_Circulant(Y,D1,D1T,D2,D2T,tau1,tau2,Tmax,maskMtx, dMtx,M,N,eps)


% Y is the n times T dimensional matrix of measurements with each column
% representing a snapshot of the acoustic wavefield.

% X is the n times d dimensional dictionary matrix. It is composed of two
% dictionaries, one corresponding to the anomalous component and one for
% the smooth component of the decomposition.

% tau1 and tau2 are the regularization constants used for
% promoting sparsity in the coefficient matrices Beta1 and Beta2,
% respectively. The coefficient matrices correspon to the first and second
% components of the demixing model.

% Tmax is the maximum number of iterations the algorithm takes.

% k is the number of atoms of the dictionary corresponding to the anomalous
% term.

tic;
T = size(Y,2);

f = zeros(Tmax,1);


XTY = [D1T(Y,M,N); D2T(dMtx, maskMtx, Y, M, N)];
% XTY = Dt(dMtx, maskMtx, Y, M, N);
k = sum(sum(maskMtx));

Alpha = zeros(size(XTY,1),T);
Beta = Alpha;
Alpha1 = Alpha(1:M*N,:);
Alpha2 = Alpha(M*N+1:end,:);

eta = 2;
Lk = 0.001;

t = 1;


for itr = 1 :Tmax
   
    XAlpha = D1(Alpha1,M,N)+D2(dMtx, maskMtx, Alpha2, M, N);
%     XAlpha = D(dMtx, maskMtx, Alpha, M, N);

    XTXAlpha = [D1T(XAlpha,M,N);D2T(dMtx, maskMtx, XAlpha, M, N)];
%     XTXAlpha = Dt(dMtx, maskMtx, XAlpha, M, N);
    
    grad = (XTXAlpha - XTY);
    
    L  = Lk ;
    flag = 1; 
    
    while flag == 1
        
            A = Alpha - (1/L)*grad;
            A1 = A(1:M*N,:);
            A2 = A(M*N+1:end,:);
    
            Beta1 = soft(A1,tau1/L);
            Beta2 = soft(A2,tau2/L);
            pLyk =  [Beta1;Beta2];
            
            XBeta = D1(Beta1,M,N) + D2(dMtx, maskMtx, Beta2, M, N);
%             XBeta = D(dMtx, maskMtx, pLyk, M, N);
            F = (1/2)*(norm(Y-XBeta,'fro'))^2 ;
            
            XAlpha = D1(Alpha1,M,N) + D2(dMtx, maskMtx, Alpha2, M, N);
%             XBeta = D(dMtx, maskMtx, Beta, M, N);
            Q = (1/2)*(norm(Y-XAlpha,'fro'))^2 + (pLyk - Alpha)'*grad ...
              + (L/2)*norm(pLyk - Alpha,'fro')^2 ;
            
            if F > Q 
               L = eta*L; 
            else
                flag = 0;
            end
            
            
    end
    
    Lk = L;
    
    A = Alpha - (1/Lk)* (XTXAlpha - XTY);
    A1 = A(1:M*N,:);
    A2 = A(M*N+1:end,:);
    
    Beta1 = soft(A1,tau1/Lk);
    p1 = tau1*sum(sum(abs(Beta1)));
    
    Beta2 = soft(A2,tau2/Lk);
    p2 = tau2*sum(sum(abs(Beta2)));
    
    Beta_0 = Beta;
    Beta = [Beta1;Beta2];
    
    t0 = t;
    t = (1+sqrt(1+4*t0^2))/2;
    
    Alpha = Beta + ((t0-1)/t) * (Beta - Beta_0);
    Alpha1 = Alpha(1:M*N,:);
    Alpha2 = Alpha(M*N+1:end,:);
    
    %if (norm(Beta,'fro') == 0)
    %   break;
        
    if  norm(Beta-Beta_0,'fro') <= eps*norm(Beta,'fro')
       break;
    end

    XBeta = D1(Beta1,M,N) + D2(dMtx, maskMtx, Beta2, M, N);
%     XBeta = D(dMtx, maskMtx, Beta, M, N);
    f(itr) = (1/2)*(norm(Y-XBeta,'fro'))^2 + p1 + p2;
%     figure(1);
%     plot(log(f(1:itr)));
%     drawnow;

end
% itr
f = f(1:itr-1);
time = toc;
end