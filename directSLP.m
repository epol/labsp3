function [Lambda,Y] = directSLP(q,L,N)
    % Transform the potential function from [ -L , + L ] to [ 0 , 2pi]
    v = @(xi) (L/pi)^2 * feval(q,(L/pi)*(xi - pi )) ;
    % Calculate the value of v on the nodes
    V = feval(v, linspace (0, 2*pi, N+1) ) ;
    % Compute the matrix for the second derivative
    D = directSLP_inner2(N) ; 
    % Resolve the problem and get the eigenvalues and eigenvectors (coordinate of the eigenfunctions in the FPM base)
    [ E, y] = directSLP_inner1(D,V) ;
    % Get the eigenfunctions as functions
    %yf = _directSLP_FPMbase(N) * y ;
    %yf = @(t) 0 ; %FIXME
    % Take back the eigenvalues of the rescaled problem to the original problem
    Lambda = (pi/L)^2 * E ; 
    %Y = @(x) feval(yf, pi/L*(x+L)) ; 
    Y = y;
end

