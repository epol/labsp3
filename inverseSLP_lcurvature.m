function alpha0 = inverseSLP_lcurvature(W,T,Beta,Sigma)
    % returns the negative of the curvature of the L-curve
    % kappa(alpha) with alpha a row
    
    % parameters
    M = length (Sigma) ;
    
    % auxiliary variables
    Beta2 = Beta .^2 ;
    Beta4 = Beta2 .^2 ;
    Beta6 = Beta2 .* Beta4 ; 
    Sigma2 = Sigma .^2 ;
    Sigma2Beta4 = Sigma2 .* Beta4 ;
    WT = W' * T ; 
    WT2 = WT.^2 ;
    
    % auxiliary functions
    % f = norm(A * Deltav - T) ^2
    f = @(alpha)  alpha.^2 .* sum ( (( (WT .* Beta2) * ones(1,length (alpha)) ) ./ (Sigma2 * ones(1,length (alpha)) + Beta2 * alpha ) ).^2 ) ;
    df = @(alpha) -2 * alpha .* sum ( ((WT2 .* Sigma2Beta4) * ones(1,length (alpha))) ./ ( Sigma2* ones(1,length (alpha)) + Beta2 * alpha ).^3 ) ;
    ddf = @(alpha) -2 * sum ( ((WT2 .* Sigma2Beta4) * ones(1,length (alpha))) .* ( Sigma2* ones(1,length (alpha)) -2* Beta2 * alpha  ) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha ).^4 ) ;
    % h = norm(D * Deltav) ^2
    h = @(alpha) sum (( ((WT .* Beta .*Sigma) * ones(1,length (alpha)) ) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha) ).^2) ;
    dh = @(alpha) -2 * sum ( ((WT2 .* Sigma2Beta4 ) * ones(1,length (alpha))) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 *alpha ).^3 ) ;
    ddh = @(alpha ) 6 * sum ( (( WT2 .* Sigma2 .* Beta6 ) * ones(1,length (alpha)) )./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha).^4 ) ;
    
    drho = @(alpha) .5 * df (alpha) ./ sqrt(f(alpha)) ;
    ddrho = @(alpha) .25 * ( 2* ddf(alpha) .* f(alpha) - (df(alpha)).^2 ) ./ ( (f(alpha)).^1.5 ) ;
    deta = @(alpha) .5 * dh (alpha) ./ sqrt(h(alpha)) ;
    ddeta = @(alpha) .25 * ( 2* ddh(alpha) .* h(alpha) - (dh(alpha)).^2 ) ./ ( (h(alpha)).^1.5 ) ;
        
    
    % (negative) curvature
    
    kappa = @(alpha) - ( drho(alpha) .* ddeta(alpha) - ddrho(alpha) .* deta(alpha) ) ./ (( (drho(alpha)).^2 + (deta(alpha)).^2 ) .^ 1.5) ;
    
    % OLD
    %kappa = @(alpha) - 2 * f(alpha) .* h(alpha) .* ( df(alpha) .* h(alpha) .*(ddh(alpha) - (dh(alpha)).^2) - f(alpha) .* dh(alpha) .* (ddf(alpha) - (df(alpha)).^2)  ) ./ ( ((df(alpha)).^2 + (dh(alpha)).^2 ) .^(1.5) ) ;
    
    
    alpha0 = fminbnd(kappa,1e-15,1e-1) 
    
    %DEBUG
    %clf
    %hold on
    alpha = linspace(1e-16, 1e1, 1000);
    %semilogy(sqrt(f(alpha)),sqrt(h(alpha)),'b')
    appa = kappa(alpha) ;
    appa = max(abs(sqrt(h(alpha))))/max(abs(appa)) * appa ;
    %plot(sqrt(f(alpha)),appa,'g')
    %DEBUG
    %semilogy(sqrt(f(alpha0)),sqrt(h(alpha0)),'r*')
    %pause
    
end
