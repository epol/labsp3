function D = inverseSLP_hermitematrix(xi)
    % compute the penalization matrix computing the second derivative on
    % the grid points after an interpolation with the hermite base
    
    % let xi be a column
    xi = xi(:) ;
    n = length(xi) ;
    
    % calculate the polynomials with the formula H_n = 2x H_(n-1) - 2(n-1) H_(n-2)
    % note: I don't use the "standard" matlab form for polynomial
    % coefficients, the i-th term of the vector is the coefficient for 
    % x^(i-1)
    % TODO: is there a way to avoid "for"? I don't think so
    hermitecoeff = zeros(n,n) ;
    hermitecoeff(1,1) = 1 ;
    hermitecoeff(2,2) = 2 ;
    for i = (3:n)
        % 2x H_(i-1)
        hermitecoeff(i,2:i) = 2*hermitecoeff(i-1,1:i-1) ;
        % -2(i-2) H_(n-2)   % this is a DAMN off-by-one becuase of matlab indexint
        hermitecoeff(i,1:i-2) = hermitecoeff(i,1:i-2) - 2*(i-2) * hermitecoeff(i-2,1:i-2) ;
    end
    
    % calculate a matrix with x_j ^i at the i,j position
    powerxi = (ones(n,1) * xi' ) .^ ( (0:n-1)' * ones(1,n) ) ;
    
    % calculate the values of each polynomial on evetry grid point
    % H(i,j) = H_i(xi_j)  with off-by-one
    H = hermitecoeff * powerxi ;
    
    % calculate the matrix for the second derivatives at the grid point
    % S2 takes care of the shifts, M for the i(i-1) coefficients
    eyen = eye(n) ; % because matlab refuses to run eye(n)(:,1:n-2)
    S2(:,3:n) = eyen (:,1:n-2) ;
    M = diag((0:n-1)) .* diag( cat( 2,[0],(0:n-2))) ;
    
    D = H * S2 * M * inv(H) ;
    
end
    
