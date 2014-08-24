function q=inverseSLP(L,Lambda,Kmax,tol)
    M = length(Lambda) ;
    N = 2*M +1 ;
    % force the Lambda vector to be a column
    Lambda = Lambda(:) ;
    % transform the eigenvalues 
    E = (L/pi)^2 * Lambda ;
    % nodes
    Xi = 2*pi/(2*M+1)*(1:M) ;
    % choose the initial potential vector 
    vk = zeros(M,1) ;  % TODO: there are a better choice?
    % we don't need to compute the differentiation every time
    D = directSLP_inner2(N-1);
    
    % main cycle
    k = 0 ;
    Deltavk = 2*vk*tol ;
    while ((norm(Deltavk) >= norm(tol * vk)) && (k <= Kmax )) %FIXME: use norms
        % extend v with symmetry
        vext(1) = 0 ;
        vext(2:M+1) = vk ;
        vext(M+1:2*M) = flipud(vk) ;
        vext(2*M+1) = 0;
        % resolve the direct problem in this case
        [ Ek, Yk] = directSLP_inner1(D,vext) ;
        Tk = Ek(1:M) - E ;
        % calculate the Jacobian matrix using the formula a_{mn} = 2(y_{n;m})^2 
        Ak = 2 * ((Yk(1:M,1:M))').^2 ;  %TODO: is this correct? I don't think so
        % TODO: approximate Deltavk with TP regularization
        % DEBUG
        condAk = cond(Ak)
        
        % SVD of Ak
        [ W , Sigma, U ] = svd(Ak) ;
        Sigmav = diag(Sigma) ;
        
        % Tikhonov regularization
        Beta = ones(M,1) ;
        kappa = inverseSLP_lcurvature(W,Tk,Beta,Sigmav) ; 
        alpha = fminbnd(kappa,1e-10,1e-1)  %TODO: something better?
        
        % Find Delta
        %Deltavk = U * (Sigmav ./ ( Sigmav .^2 + alpha * Beta .^2) .* (W' *Tk)  ) ;
        Deltavk = ridge(Tk,Ak,Beta) ; 
        
        %DEBUG
        before = norm(Ak*vk - Tk )
        vk = vk - Deltavk ;
        %DEBUG
        after = norm(Ak*vk - Tk)
        
        vk = vk - Deltavk ;
        k = k+1
        %Deltavk >= tol*vk 
    end
    q = vk ;
end
