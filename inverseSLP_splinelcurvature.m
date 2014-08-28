function alpha0 = inverseSLP_splinelcurvature(W,T,Beta,Sigma)
    M = length (Sigma) ;
    
    % auxiliary variables
    Beta2 = Beta .^2 ;
    Beta4 = Beta2 .^2 ;
    Beta6 = Beta2 .* Beta4 ; 
    Sigma2 = Sigma .^2 ;
    Sigma2Beta4 = Sigma2 .* Beta4 ;
    WT = W' * T ; 
    WT2 = WT.^2 ;

    % f = norm(A * Deltav - T) ^2
    f = @(alpha)  alpha.^2 .* sum ( (( (WT .* Beta2) * ones(1,length (alpha)) ) ./ (Sigma2 * ones(1,length (alpha)) + Beta2 * alpha ) ).^2 ) ;
    % h = norm(D * Deltav) ^2
    h = @(alpha) sum (( ((WT .* Beta .*Sigma) * ones(1,length (alpha)) ) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha) ).^2) ;

    % interpolation points
    alpha = logspace(-16,0,50) ; 

    x = sqrt(f(alpha)) ;
    y = sqrt(h(alpha)) ;
    
    % blur the interpolation points 
    xx = .2 * ( x(1:end-2) + 3*x(2:end-1) + x(3:end)) ;
    yy = .2 * ( y(1:end-2) + 3*y(2:end-1) + y(3:end)) ;
    pp = spline (xx,yy) ;
    
    % derivates of the spline
    dpp = pp ; 
    dpp.coefs = dpp.coefs *  diag([ 3 2 1 ] , 1) (:,2:end) ;
    dpp.order = dpp.order -1 ;
    ddpp = dpp ; 
    ddpp.coefs = ddpp.coefs *  diag([ 2 1 ] , 1) (:,2:end) ;
    ddpp.order = ddpp.order -1 ;
    
    % (negative) curvature of the spline
    kappa = @(x) - ppval(ddpp,x) ./ (1 + ppval(dpp,x).^2).^(1.5) ;
    
    % a (limited) domain for pp
    x0 = sqrt(f(1e-15)) ;
    x1 = sqrt(f(1e-1));
    
    % find the corner
    xmin = fminbnd(kappa,x0,x1) ;

    % take back the coordinate of the corner to the orginal alpha parameter
    funxzero = @(alpha) sqrt(f(alpha)) - xmin ; 
    alpha0 = fzero(funxzero, [1e-15,1e-1]) ;
    
    
    % DEBUG
    %minimovettoredebug = [x0,xmin,x1]
    clf
    hold on
    plot(sqrt(f(alpha)),sqrt(h(alpha))) 
    plot(x,ppval(pp,x),'g')
    plot(xmin,ppval(pp,xmin),'r*') 
    pause
    


end
