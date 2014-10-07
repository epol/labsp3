function q=inverseSLP(L,Lambda,Kmax,tol,v0)

    %NEW: regtools
    addpath('./regu')
    
    if nargin < 4
        error('Too few parameters')
    end
    if nargin < 5
        v0 = zeros(length(Lambda),1)
    end
    
    
    M = length(Lambda) ;
    N = 2*M +1 ;
    % force the Lambda vector to be a column
    Lambda = Lambda(:) ;
    % transform the eigenvalues 
    E = (L/pi)^2 * Lambda ;
    % nodes
    Xi = 2*pi/(2*M+1)*(1:M) ;
    Xi = Xi(:) ;
    % choose the initial potential vector 
    vk = v0 ;  % TODO: there are a better choice? They suggest a q + b, so maybe this should be an argument for this function
    % we don't need to compute the differentiation every time
    D = directSLP_inner2(N-1);
    
    %compute a matrix that extends the potential with symmetry
    I = eye(M) ;
    extender = [ I(1:end,:) ; zeros(1,M) ; flipud(I(1:end,:)) ];
    
    size(D)
    size(extender)
    
    % main cycle
    k = 0 ;
    Deltavk = 2*vk*tol ; % a way to implement a "do while"
    while ((norm(Deltavk) >= norm(tol * vk)) && (k <= Kmax )) %FIXME: additional stopping condition if we know noise level
        % extend v with symmetry
        vext(1) = 0 ;
        vext(2:M+1) = vk ;
        vext(M+1:2*M) = flipud(vk) ;
        vext(2*M+1) = 0;
        % resolve the direct problem in this case
        [ Ek, Yk] = directSLP_inner1(D,vext) ;
        
        Tk = Ek(1:M) - E ;
        
        %DEBUG
        normTk = norm(Tk)
        
        % calculate the Jacobian matrix using the formula a_{mn} = 2(y_{n;m})^2 
        Ak = 2 * ((Yk(1:M,1:M))').^2 ;
        
        % SVD of Ak
        %[ W , Sigma, U ] = svd(Ak) ;
        %Sigmav = diag(Sigma) ;
        
        % Tikhonov regularization
        %Beta = ones(M,1) ; %TODO: D = W * Beta * U' should be the second order hermite differntation matrix to favor smooth solutions
        %Beta = min(100*ones(M,1),ones(M,1)./Sigmav) 
        
        % calculate the regularization parameter with the L-curve method
        %alpha = inverseSLP_lcurvature(W,Tk,Beta,Sigmav) 
        
        % OLD
        %alpha = fminbnd(kappa,1e-7,1e1)  %TODO: something better?
        
        % NEW try with library
        %% computer the differntation matrix
        %%% hermite
        %scale = herroots(M);
        %scale = scale(M)/pi;
        %[r,DiffMat] = herdif(M,2,scale) ;
        %%% poly
        %DiffMat = poldif(Xi,3);
        %DiffMat = DiffMat(:,:,3) ;
        %%% orginal differntiation matrix from the original problem
        DiffMat = D * extender;
        DiffMat = DiffMat(1:M);
        
        %% svd and gsvd 
        %[ W, Sigma, U ] = csvd(Ak);
        size(Ak)
        size(DiffMat)
        [ WW, SigmaM, XX] = cgsvd(Ak,DiffMat) ;
        %% find the optimal parameter
        [reg_corner,rho,eta,reg_param] = l_curve(WW,SigmaM,Tk,'Tikh') ;
        %% compute the transformation
        Deltavk = tikhonov(WW,SigmaM,XX,Tk,reg_corner) ;
        pause
        
        
        % Find Delta
        %Deltavk = U * (Sigmav ./ ( Sigmav .^2 + alpha * Beta .^2) .* (W' *Tk)  ) ;
        
        % Calculate the new point        
        vk = vk - Deltavk ;
        k = k+1 ;
        %Deltavk >= tol*vk 
    end
    
    % TODO give a nice rescaled output, or a function.
    q = vk ;
    
end

