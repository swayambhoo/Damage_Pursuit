function [ X1, X2, Obj, time ] = Damage_Pursuit_v4( Data , varargin)

% Damage Pursuit solves optimization problems of the form 
%
%              minimize Loss(y;x1,x2)+f1(x1)+f2(x2)        (P1)
%
% with respect to x1 and x2, where y is the observation vector/matrix.
% Assume y is linearly decomposable as the sum of two linear functions of
% x1 and x2, i.e.
%
%                       y = D1 * x1 + D2 * x2, 
%
% with D1 and D2 being predefined dictionaries.Then, an appropriate choice 
% for the loss function is the residual sum of squares 
%
%            Loss(y;x1,x2) = 1/2 * ||y - (D1 * x1 + D2 * x2)||^2.
%
% The variables x1 and x2 are assumed to have sparse structures imposed by 
% penalty functions f1 and f2.

% Damage pursuit is based on the l1-norm penalty
% 
%                    f1(x1) = tau1 * ||x1||_1 and 
%
% for x1 and l1/l2-norm penalty for x2, i.e.
%
%                  f2(x2) = tau2 * sum_i(||x2,i||_2)  
%
% where indexing of x2 depends on how the spatial/temporal grouping of the 
% pixels is performed. The groups are assumed to be non-overlapping here. 

% -------------------------------------------------------------------------

%   ============  Required Inputs  =========== 

% Data: Data cube whose slices store images acquired
% at different time instants.

%  ============  Optional Inputs   ===========
% 
%   D1_type : must be one of {0,1,2}
%       
%        0  :  D1 is set to be a function handle which computes the two
%              dimensional inverse Discrete Cosine Transforms (IDCT) of the
%              frames and stores the results in a matrix whose columns 
%              correspond to different time frames.
%       

%        1  :  D1 is set to be a function handle which computes the two
%              dimensional inverse Fast Fourier Transforms (IFFT) of the
%              frames and stores the results in a matrix whose columns 
%              correspond to different time frames.
%
%        2  :  Function_handle will be provided by the user. If D1 is
%              provided by the user, D1T, the function handle to compute the 
%              inverse transform of D1 should also be provided.             
%
%   Default: 0, i.e. DCT transform will be used by default.
%   
%   D1: A dictionary for the smooth component or a function handle which
%   computes products of the form D1 * x1.
%
%   D1T: If D1 is a function handle, then D1T should also be inputed as a
%   handle to a function that computes products of the form D1T * x.
%   If D1 is a matrix, D1T does not need to be inputed.

%   SNR : Signal to Noise ratio. The added noise will be i.i.d. Gaussian.
%   Default: SNR = inf. No niose is added to the measurements.

%  'tau1' : Regularization parameter in f1.
%           Default: tau1 = Sig * sqrt(log(Tmp_Rsl * dim1 * dim2 / delta1 ));

%  'tau2' and 'tau3' : Regularization parameters in f2. tau3 is multiplied with 
%   the l1/l2 norm of the marginal pixels of the image, whereas tau2 is used 
%   for the rest of the pixels.
%   Default:  
%   Default: tau3 = 2 * tau2.

%  'Dws_fac' : Downsampling factor.
%   Default: Dws_fac = 1.

%  'Spt_Rsl' : Spatial resolution is used to group the images' pixels. 
%   It gives the side length of the pixel cubes covering images. 
%   Default: Spt_Rsl = 1, i.e. pixels are not grouped.
 
%  'Tmp_Rsl' : Temporal resolution determines the number of image frames  
%   for which the optimization problem is run simultaneously. 
%   Default: Tmp_Rsl = 1, i.e. the frames are solved separately.

%  'Margine' : Number of marginal cubes of pixels whose l1-l2 norms are 
%   multiplied with tau3. 
%   Default: Margine = 0.

%  'eps' : Required accuracy for solving the optimization problem P1.
%   Default: eps = 1e-4.

%  'Initial' : must be one of {0,1,2}
%               0 : The problem of the first frame will be initialized at 
%                   zero and the while the next ones will be warm started.
%
%               1 : Zero Inizialization for all the frames.
%
%               2 :  The problem of the first frame will be initialized
%               randomly but the next ones will be warm started.
%
%               3 : Random Initialization for all the frames.
%
%               Default = 0;

%  'Max_Itr' : Maximum number of iterations of the algorithm.
%   Default: Max_Itr = 500.

%  'Trm_Crt' : Termination criterion to use:
%              0 : norm of the distance of two consecutive iterates relative
%               to the norm of one iterate is used as the termination
%               criterion.
%              1 : the norms of the residuals are used as for the stopping
%               cretirion.
%   Default: Trm_Crt = 0.

% 'Plot_flag' : If this flag is one, the reconstructed coefficients will be
%  plotted.
%  Default: Plot_flag = 0.
 
%  ============  Outputs   =================

%  =========================================

T = size(Data,3);
dim1 = size(Data,1);
dim2 = size(Data,2);


for t = 1 : T
    y(:,t) = reshape(Data(:,:,t),dim1*dim2,1);
end

%   Set default values of variables.
Dws_fac = 1;
Spt_Rsl = 1;
Tmp_Rsl = 1;
Margine = 0;
eps = 1e-6;
Trm_Crt = 0;
Max_Itr = 100;
Initial = 0;
Plot_flag = 0;
D1_type = 0;
D2_type = 0;
thr = 1e-4;
sgm = 1e-5;
% test for number of required inputs.

if (nargin-length(varargin)) ~= 1
    error('Please input the data cube to the software');
end

if (rem(length(varargin),2)==1)
    error('Optional inputs should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'TAU1'
                tau1 = varargin{i+1};
            case 'TAU2'
                tau2 = varargin{i+1};
            case 'TAU3'
                tau3 = varargin{i+1};
            case 'DWS_FAC'
                Dws_fac = varargin{i+1};
            case  'SPT_RSL'
                Spt_Rsl = varargin{i+1};
            case 'TMP_RSL'
                Tmp_Rsl = varargin{i+1};
            case 'MARGINE'
                Margine = varargin{i+1};
            case 'EPS'
                eps = varargin{i+1};
            case 'TRM_CRT'
                Trm_Crt = varargin{i+1};
            case 'D1_TYPE'
                D1_type = varargin{i+1};
            case 'D1'
                D1 = varargin{i+1};
            case 'D1T'
                D1T = varargin{i+1};
            case 'D2_TYPE'
                D2_type = varargin{i+1};
            case 'D2'    
                D2 = varargin{i+1};
            case 'D2T'
                D2T = varargin{i+1};
            case 'THR'
                thr = varargin{i+1};
            case 'SGM'
                sgm = varargin{i+1};
            case 'LX'
                Lx = varargin{i+1};
            case 'LY'
                Ly = varargin{i+1};
            case 'MAX_ITR'
                Max_Itr = varargin{i+1};
            case 'INITIAL'
                Initial = varargin{i+1};
            case 'PLOT_FLAG'
                Plot_flag = varargin{i+1};
                coeffsmat2 = 0;
            otherwise
                % The parameter string is wrong.
                error(['Unrecognized Input: ''' varargin{i} '''']);
        end
    end
end

%%%%%%%%%%%%%%%%%%
% First normalize the matrix of measurements.
% w = y./max(max(abs(y)));
w = y;
% s1 and s2 give the image size and the number of time frames, respectively.
s1 = size(w,1);
s2 = size(w,2);

% The images are downsampled at the specified rate Dws_fac. 
if Dws_fac ~= 1
    
    wnew = zeros(0,0);
    
    for ii = 1:s2
        
        wtemp = w(:,ii);
        wtemprs = reshape(wtemp, dim1, dim2);
        wds = wtemprs(1:Dws_fac:end, 1:Dws_fac:end);
        dim1new = size(wds,1);
        dim2new = size(wds,2);
        wnew(:,ii) = wds(:);
        
    end
    
    w = wnew;
    s1 = size(w,1);
    s2 = size(w,2);
    
    dim1 = dim1new;
    dim2 = dim2new;
end

% The transformations are specified here.

switch D1_type
    
    case 0
        
        D1 = @(x,dim1,dim2) Dmult(x,dim1,dim2);
        D1T = @(x,dim1,dim2) DTmult(x,dim1,dim2);
        
    case 1
        
        D1 = @(x,dim1,dim2) Fmult(x,dim1,dim2);
        D1T = @(x,dim1,dim2) FTmult(x,dim1,dim2);
        
    case 2
        
        if isa(D1,'function_handle')
            
            try
                D1(randn(size(y)),dim1,dim2);
            catch
                error('The function handle D1 is not working');
            end
            
            if ~exist('D1T')
                error('If D1 is a function handle, D1T should also be entered');
            else   
                try
                    D1T(randn(size(y)),dim1,dim2);
                catch
                    error('The function handle D1T is not working');
                end
                
            end
            
        else
            
            D1T = @(x,dim1,dim2) D1' * x;
            D1 = @(x,dim1,dim2) D1 * x;
            
            
        end
end

% The second transformation is specified here.

switch D2_type
    
    case 0
        
        % In this case the dictionary for the second component will be the
        % identity dictionary. We will use the alternating algoright in
        % this case in order to solve the corresponding sub-problems.
        
        D2T = @(x) x;
        D2 = @(x) x;
        
        
        % The pixels in each image form groups of the size specified by Spt_Rsl.
        
        n1 = floor(dim1/Spt_Rsl);
        n2 = floor(dim2/Spt_Rsl);
        
        if Spt_Rsl ~= 1
            
            % Dimensions are rounded to multiples of Spt_Rsl.
            
            dim1new = Spt_Rsl * n1;
            dim2new = Spt_Rsl * n2;
            
            if dim1 ~= dim1new || dim2 ~= dim2new
                
                for j = 1 : s2
                    
                    wtemp = reshape(w(:,j),dim1,dim2);
                    wtemp2 = wtemp(1:dim1new,1:dim2new);
                    wnew(:,j) = wtemp2(:);
                end
                
                w = wnew;
                dim1 = dim1new;
                dim2 = dim2new;
            end
            
        end
        
        % Here we compute a matrix called Index which specifies how pixels are
        % indexed in groups.
        
        tmp1 = reshape(1 : (n1 * n2) , n1 , n2);
        tmp2 = ones(Spt_Rsl,Spt_Rsl);
        tmp3 = kron(tmp1 , tmp2);
        [~ , Index] = sort(tmp3(:));
        
        % Tau is a matrix whose entries give the coefficients for the penalty
        % function f2.
        
        A1 = kron([1:n1]',ones(1,n2));
        A2 = kron(1:n2,ones(n1,1));
        
        if ~exist('tau3','var')
            tau3 = 2 * tau2;
        end

        Tau = tau2 * ones(n1 , n2) + (tau3-tau2) * max( ...
            (min (A1 , n1 - A1 + 1) <= Margine) ,...
            (min (A2 , n2 - A2 + 1) <= Margine));
        
        Tau = Tau(:);
        
    case 1
        
        % In this case the columns of the dictionary are set to vectorized
        % versions of the 2D Marr wavelets centered at interior points of
        % the structure.
       
        for i = 1 : length(sgm)
            
            [d,mask] = MarrWvlt_Dct(Lx,Ly,dim1,dim2,sgm(i),thr(i));
            dMtx(:,i) = d;
            maskMtx(:,i) = mask;
            
        end
        
        D2  = @(dMtx, maskMtx, X, dim1, dim2) Cmult...
            (dMtx, maskMtx, X, dim1, dim2);
        D2T = @(dMtx, maskMtx, Y, dim1, dim2) CTmult...
            (dMtx, maskMtx, Y, dim1, dim2);
        
    case 2
        
        if isa(D2,'function_handle')
            
            try
                D2(randn(size(y)),dim1,dim2);
            catch
                error('The function handle D2 is not working');
            end
            
            if ~exist('D2T')
                error('If D2 is a function handle, D2T should also be entered');
            else   
                try
                    D2T(randn(size(y)),dim1,dim2);
                catch
                    error('The function handle D2T is not working');
                end
                
            end
            
        else
            
            D2T = @(x,dim1,dim2) D2' * x;
            D2 = @(x,dim1,dim2) D2 * x;
            
        end
end

%%%%%%%%%%%%%%%%





% X1_0 and X2_0 are the initial points for the algorithm.

if Initial <= 1
    X1_0 = zeros(dim1 * dim2 , Tmp_Rsl);
    X2_0 = zeros(dim1 * dim2 , Tmp_Rsl);
else
    X1_0 = randn(dim1 * dim2 , Tmp_Rsl);
    X2_0 = randn(dim1 * dim2 , Tmp_Rsl);
end

% ii is the iteration counter.
ii = 1;

while ii <= ( s2 - Tmp_Rsl + 1 )
    
    % Get current frames
    yy = w(: , ii : ii + Tmp_Rsl - 1); 
    
    if D2_type == 0
        
    [XTmp1, XTmp2 , Obj , time] = Alternating_Algorithm ...
        (yy , D1 , D1T , dim1 , dim2 , Index , X1_0 , X2_0 , ...
         tau1 , Tau , eps , Spt_Rsl , Max_Itr , Trm_Crt);
    
    elseif D2_type == 1
        
       [XTmp1, XTmp2, Obj , time] = FISTA_Circulant...
        (yy, D1,D1T, D2, D2T, tau1, tau2, Max_Itr, maskMtx, ...
         dMtx, dim1, dim2, eps);
   
    end
    
    X1(: , ii : ii + Tmp_Rsl - 1) = XTmp1;
    X2(: , ii : ii + Tmp_Rsl - 1) = XTmp2;
    
     if Plot_flag == 1
         
         x1 = D1(XTmp1,dim1,dim2);
         % Identify the support of the coefficients 
         coeffsmat = double(abs(XTmp2)>5*tau2);

         for j = 1 : Tmp_Rsl
             
             figure(2);
             subplot(2,3,1);
             imagesc( reshape(yy(:,j),dim1,dim2) ); colormap jet
%              set(gca,'YDir','normal');
             axis equal, axis off
             title('Original mixture')
             
             subplot(2,3,2);
             imagesc( reshape(x1(:,j),dim1,dim2) ); colormap jet
%              set(gca,'YDir','normal');
             axis equal, axis off
             title('Recovered Smooth')
             
             subplot(2,3,3);
             imagesc( abs(reshape(XTmp2(:,j),dim1,dim2)) ); colormap jet
%              set(gca,'YDir','normal');
             axis equal, axis off
             title('Recovered Sparse')
             
             subplot(2,3,5);
             imagesc( reshape(abs(coeffsmat(:,j)),dim1,dim2) ); colormap jet
%              set(gca,'YDir','normal');
             axis equal, axis off
             title('Recovered coefficients (support)')
             
             coeffsmat2 = coeffsmat2 + coeffsmat;
             subplot(2,3,6);
             imagesc( reshape(abs(coeffsmat2(:,j)),dim1,dim2) ); 
%              set(gca,'YDir','normal'); 
             colormap jet
             axis equal, axis off             
             title('Sparse Persistence Map');
             
             drawnow;
             
         end
     end
     
     switch Initial
         
         case 0 , 2;
             X1_0 = XTmp1;
             X2_0 = XTmp2;
         case 1
             X1_0 = zeros(dim1 * dim2 , Tmp_Rsl);
             X2_0 = zeros(dim1 * dim2 , Tmp_Rsl);
         case 3
             X1_0 = randn(dim1 * dim2 , Tmp_Rsl);
             X2_0 = randn(dim1 * dim2 , Tmp_Rsl);
             
     end
           
     ii = ii + Tmp_Rsl

end

end

