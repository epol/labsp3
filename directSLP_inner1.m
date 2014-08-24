function [ E, Y] = _directSLP_inner1(D,V)
    % D is D^(2).
    % v is the rescaled potential function
    N = length(V) - 1;
    % xi = linspace (0, 2*pi, N+1);
    % xir = xi(2:end-1);   % we are not interested in the values on 0 and 2pi
    vir = V(2:end-1) ; % feval(v,xir);
    
    M = - D(2:end-1,2:end-1) + diag(vir);
    [Y,D] = eig(M);
    E = diag(D) ;
    E = E(:) ;
end
