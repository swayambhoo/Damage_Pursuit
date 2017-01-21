clc; clear all; close all;

M = 32;
N = 32;

Lx = 0.75;
Ly = 0.75;

sigma = 1e-2;
thr = 1e-5;

for i = 1 : length(sigma)
    
    [d,mask] = MarrWvlt_Dct(Lx,Ly,M,N,sigma(i),thr(i));
    dMtx(:,i) = d;
    maskMtx(:,i) = mask;
    
end

% col = d;
% row = [col(1) fliplr(col(2:end))];
% D = toeplitz(col,row);
% 
% for i = 1 : size(D,2)
%     if mask(i) == 1
%         figure(10);
%         imagesc(reshape(D(:,i),M,N));
%         drawnow;
%     end
% end

n = sum(maskMtx);

% Here we generate a toy example problem to test Damage Pursuit.
T = 1;
% X1 is the coefficient for the smooth term.
X1 = zeros(M*N,T);
X1([1 20 40],:) = 1;

% X2 is the coefficient for the anomaly term.
X2 = zeros(n,T);
X2([10 100 200],:) = 1;

% X is the coefficient matrix.
X = [X1 ; X2];

D1_Type = 1;

if D1_Type == 0;
    Y = Cmult(dMtx, maskMtx, X2, M, N) + Dmult(X1,M,N);
elseif D1_Type == 1;
    Y = Cmult(dMtx, maskMtx, X2, M, N) + (Fmult(X1,M,N));
end

tau1 = .1;
tau2 = .1;
Tmax = 10000;


[X1hat , X2hat, Obj , time ] = Damage_Pursuit_v2...
    (Y,M,N,'sgm',sigma,'thr',thr,'Lx',Lx,'Ly',Ly,'D2_Type',1,...
     'D1_Type',D1_Type,'tau1',tau1,'tau2',tau2);

figure(1);
plot(X1);
hold on;
plot(real(X1hat),'red');

figure(2);
plot(X2);
hold on;
plot(X2hat,'red');

if D1_Type == 0
    
    figure(3);
    Yhat = Cmult(dMtx, maskMtx, X2hat ,M, N)+Dmult(X1hat,M,N);
    imagesc(reshape(Yhat,M,N));
    
    figure(4);
    imagesc(reshape(Y,M,N));
    
end

