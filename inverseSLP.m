function q=inverseSLP(L,Lambda,Kmax,tol,v0)

    %NEW: regtools
    addpath('./regu')
    addpath('./dmsuite')
    
    if nargin < 4
        error('Too few parameters')
    end
    if nargin < 5
        v0 = zeros(length(Lambda),1);
    end
    v0 = v0(:);
    
    
    M = length(Lambda) ;
    N = 2*M +1;
    % force the Lambda vector to be a column
    Lambda = Lambda(:) ;
    % transform the eigenvalues 
    E = (L/pi)^2 * Lambda ;
    % nodes
    Xi = 2*pi/(N)*(1:M) ;
    Xi = Xi(:) ;
    % choose the initial potential vector 
    vk = v0 ;  % TODO: there are a better choice? They suggest a q + b, so maybe this should be an argument for this function
    % we don't need to compute the differentiation every time
    D = directSLP_inner2(N);
    
    %compute a matrix that extends the potential with symmetry
    extender = [ zeros(1,M) ; eye(M) ;  flipud(eye(M)) ; zeros(1,M)];
    
    %size(D)
    %size(extender)
    
    % NEW try with library
    %% computer the differntation matrix
    %%% hermite
    %scale = herroots(M);
    %scale = scale(M)/pi;
    %scale=1
    %[r,PenMat] = herdif(N+1,2,scale) ;
    %%% poly
    %PenMat = poldif(Xi,2);
    %roots = herroots(M);
    %roots = pi / roots(M) * roots;
    %PenMat = poldif(herroots(M),2);
    %PenMat = PenMat(:,:,2) ;
    %PenMat = PenMat *extender;
    %PenMat = PenMat(2:M+1,:)
    %%% orginal differntiation matrix from the original problem
    %PenMat = D * extender;
    %PenMat = PenMat(2:M+1,:);
    %PenMat = D(2:M+1,2:M+1);
    %PenMat = eye(M);
    PenMat = hermite_differentation_matrix(M);
    
    
    % main cycle
    k = 0 ;
    Deltavk = 2*vk*tol ; % a way to implement a "do while"
    while ((norm(Deltavk) >= norm(tol * vk)) && (k <= Kmax )) %FIXME: additional stopping condition if we know noise level
        % extend v with symmetry
        %vext(1) = 0 ;
        %vext(2:M+1) = vk ;
        %vext(M+1:2*M) = flipud(vk) ;
        %vext(2*M+1) = 0;
        vext = extender * vk;
        % resolve the direct problem in this case
        [ Ek, Yk] = directSLP_inner1(D,vext) ;
        
        Tk = Ek(1:M) - E ;
        nT(k+1) = norm(Tk);
        
        %DEBUG
        normTk = norm(Tk)
        
        % calculate the Jacobian matrix using the formula a_{mn} = 2(y_{n;m})^2 
        Ak =  2*(((Yk(2:M+1,1:M))').^2) ;
        
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
        
        
        %% svd and gsvd 
        %[ W, Sigma, U ] = csvd(Ak);
        %sizeak = size(Ak)
        %sizePenMat = size(PenMat)
        [ WW, SigmaM, XX, VV] = cgsvd(Ak,PenMat) ;
        %% find the optimal parameter
        %[reg_corner,rho,eta,reg_param] = l_curve(WW,SigmaM,Tk,'Tikh') ;
        %reg_corner = min(reg_corner,1e0);
            
        reg_param = logspace(0,-3,600);
        [x,rho,eta] = tikhonov(WW,SigmaM,XX,Tk,reg_param,vk);
        [reg_corner,rho_c,eta_c] = l_corner (rho,eta,reg_param,WW,SigmaM,Tk,'Tikh');
        %reg_corner 
        %[reg_corner,rho,eta,reg_param] = lcurve2(WW,SigmaM,Tk);
        reg_params(k+1) = reg_corner ;
        reg_corner
        
        %% compute the transformation
        Deltavk = tikhonov(WW,SigmaM,XX,Tk,reg_corner,vk) ;
        %pause
        
        
        % Find Delta
        %Deltavk = U * (Sigmav ./ ( Sigmav .^2 + alpha * Beta .^2) .* (W' *Tk)  ) ;
        
        % Calculate the new point        
        vk = vk - Deltavk ;
        k = k+1 ;
        %Deltavk >= tol*vk 
        %plot(vk);
        %pause
    end
    clf
    hold on
    semilogy(reg_params)
    semilogy(nT/max(nT)*max(reg_params),'r')
    hold off
    q = extender(2:end-1,:) * ((pi/L)^2 * vk )  ;
    
end

