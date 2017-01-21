function [ X1 , X2 , f , time ] = Alternating_Algorithm_v2( Y , D , DT , dim1 , dim2 , I , X1 , X2 , tau1 , Tau , eps , Spt_Rsl , Max_Itr , Trm_Crt)

%   This function solves optimization problems of the following form
%
%              min  Loss(Y,D;X1,X2) + f1(X1) + f2(X2) 
%
%   with respect to X1 and X2, where the first term is a quadratic loss, 
%
%             Loss(Y,D;X1,X2) = 1/2 || Y - D * X1 - X2 ||^2
%
%   and f1 and f2 are the regularization functions defined as
%
%             f_1(X1) = tau1 * ||X1||_1
%
%             f_2(X2) = tau2 * sum_{i: Not Marginal}||X2_i||_2 + ...
%                       tau3 * sum_{  i:   Marginal}||X2_i||_2.
%
%   Minimizing the function f1 promotes the sparsity of X1, whereas having
%   f2 in the objective encourages the sparsity of X2 in the group level.
%   The matrix I determines how the grouping of the pixels is done.
%
%
%
%                 ================  Outputs   =================
%
%   'X1' : The optimal X1 solution.
% 
%   'X2' : The optimal X2 solution.
%
%   'f' : The sequence of function values at different iterations.
%
%   'time' : The time used for solving the problem.
%
%                ==============================================

tic;

% Precompute DT(y) since it will be used alot.
DTy = DT(Y,dim1,dim2);
% Compute the initialization of DT(X2).
DTx2 = DT(X2,dim1,dim2);

% First-order optimality condition with respect to X1 implies that 
% norm(g1-X1,inf) <= tau1.
g1 = DTy - DTx2;

% Primal residuals w.r.t. X1 and X2 are stored in Res1 and Res2.
Res1 = zeros(0,0);
Res2 = zeros(0,0);
Res2Tmp = zeros(size(X2,1)/(Spt_Rsl^2),1);

% f stores the function values of different iterations.
f = zeros(0,0);

% Distance of consecutive iterates of X1 and X2 are stored in
% Err1 and Err2, respectively.
Err1 = zeros(0,0);
Err2 = zeros(0,0);

% ii is the iteration counter.
ii = 1;
T = size(Y,3);

d = Spt_Rsl^2*T;
K = size(X2,1)/(Spt_Rsl^2);

[~,II] = sort(I);

while ii <= Max_Itr
    
    % Current iteration of X1 is saved as X1_0.
    X1_0 = X1;
    % X1 is simply computed by using the soft-thresholding operator.
    X1 = soft(g1 , tau1);
    
    Dx1 = D(X1,dim1,dim2);
    
    % Reg saves the values of the regularization function f_1(X1)+f_2(X2) at
    % the current iteration.
    Reg = tau1 * sum(sum(abs(X1)));
    
    % Current iteration of X2 is saved as X2_0.
    X2_0 = X2;
    
    % The following loop uses the block-coordinate descent method to
    % update the groups of pixels in X2.
    
    % First-order optimality condition with respect to X_2(Index,:), i.e. 
    % the subset of the entries of X_2 indexed by Index, says
    % that norm(g2+X_2(Index,:)) <= Tau(jj).
    g2 = Dx1 - Y;
    
    tmp1 = g2(I,:);
    tmp1 = tmp1';
    tmp1 = tmp1(:);
    tmp2 = reshape(tmp1,d,K);
    scl = sqrt(sum(conj(tmp2).*tmp2,1));
    scl = scl';
    scl2 = max(0,1-Tau./scl);
    tmp3 = (ones(d,1)*scl2') .* -tmp2;
    Reg = Reg + sum(Tau' .* sqrt(sum(conj(tmp3).*tmp3,1)));
    tmp3 = tmp3(:);
    tmp4 = reshape(tmp3,T,length(tmp3)/T);
    tmp4 = tmp4';
    X2 = tmp4(II,:);
    
    
    DTx2 = DT(X2,dim1,dim2);
    
    g1 = DTy - DTx2;
    
    % The parameters that will stop the algorithm are computed here.
    switch Trm_Crt
        
        case 0
            Err1(ii) = norm(X1 - X1_0,'fro');
            Err2(ii) = norm(X2 - X2_0,'fro')/(norm(X2,'fro')+eps);
            
            if Err1(ii) <= (eps * norm(X1,'fro')) ...
            && Err2(ii) <= (eps * norm(X2,'fro'))
                
                break;
                
            end
            
        case 1
            
            Res1(ii) = norm((g1-X1)/tau1,'inf')-1;
            Res2(ii) = norm(Res2Tmp./Tau,'inf')-1;
            
            if Res1(ii) < eps && Res2(ii) < eps
                
                break;
                
            end
            
    end

    f(ii) = (1/2) * (norm(Dx1 + X2 - Y , 'fro') ^ 2) + Reg;
    
    ii = ii + 1;
end


time = toc;
end

