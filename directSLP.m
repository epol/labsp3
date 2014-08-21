function [Lambda,Y] = directSLP(q,L,N)
    % Transform the potential function from [ -L , + L ] to [ 0 , 2pi]
    v = @(xi) (L/pi)^2 * feval(q,(L/pi)*(xi - pi )) ;
    % Compute the matrix for the second derivative
    D = _directSLP_inner2(N) ; 
    % Resolve the problem and get the eigenvalues and eigenvectors (coordinate of the eigenfunctions in the FPM base)
    [ E, y] = _directSLP_inner1(D,v) ;
    % Get the eigenfunctions as functions
    %yf = _directSLP_FPMbase(N) * y ;
    yf = @(t) 0 ; %FIXME
    % Take back the eigenvalues and eigenfunctions of the rescaled problem to the original problem
    Lambda = (pi/L)^2 * E ; 
    Y = @(x) feval(yf, pi/L*(x+L)) ; 
end
    
