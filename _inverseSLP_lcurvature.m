function kappa = _inverseSLP_lcurvature(W,T,Beta,Sigma)
    % kappa(alpha) with alpha a row
    
    %TODO: controllare TUTTO
    
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
    f = @(alpha)  alpha.^2 .* sum ( (( (WT .* Beta2) * ones(1,length (alpha)) ) ./ (Sigma2 * ones(1,length (alpha)) + Beta2 * alpha ) ).^2 ) ;
    df = @(alpha) -2 * alpha .* sum ( ((WT2 .* Sigma2Beta4) * ones(1,length (alpha))) ./ ( Sigma2* ones(1,length (alpha)) + Beta2 * alpha ).^3 ) ;
    ddf = @(alpha) -2 * sum ( ((WT2 .* Sigma2Beta4) * ones(1,length (alpha))) .* ( Sigma2* ones(1,length (alpha)) + Beta2 * alpha - 3 * ones(M,1) *alpha ) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha ).^4 ) ;
    h = @(alpha) sum (( ((WT .* Beta .*Sigma) * ones(1,length (alpha)) ) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha) ).^2) ;
    dh = @(alpha) -2 * sum ( ((WT2 .* Sigma2Beta4 ) * ones(1,length (alpha))) ./ ( Sigma2 * ones(1,length (alpha)) + Beta2 *alpha ).^3 ) ;
    ddh = @(alpha ) -2 * sum ( (( WT2 .* Sigma2 .* Beta6 ) * ones(1,length (alpha)) )./ ( Sigma2 * ones(1,length (alpha)) + Beta2 * alpha).^4 ) ;
    
    % result
    kappa = @(alpha) 2 * f(alpha) .* h(alpha) .* ( df(alpha) .* h(alpha) .*(ddh(alpha) - (dh(alpha)).^2) - f(alpha) .* dh(alpha) .* (ddf(alpha) - (df(alpha)).^2)  ) ./ ( ((df(alpha)).^2 + (dh(alpha)).^2 ) .^(1.5) ) ;
    
end
