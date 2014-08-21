function D = _directSLP_inner2(N)
    rif = ones(N+1,1)*(0:N) ; 
    rif = rif' - rif ; 
    % rif is now a (N+1)x(N+1) matrix with rif(k,j)=k-j
    if mod(N,2)
        % N even
        D = csc((pi/N)*rif) ;                   % D(k,j) = csc((k-j)*pi/N)
        D = D.^2 ;                              % D(k,j) = csc^2((k-j)*pi/N)
        D = - .5 * (-1).^rif .*D ;              % D(k,j) = -6/12 (-1)^(k-j)*csc^2((k-j)*pi/N)
        D(1:(N+2):end) = - (N^2 +2)/12 ;        % change elements on the diag
    else
        % N odd
        D = csc((pi/N)*rif) .* cot((pi/N)*rif); % D(k,j) = csc((k-j)*pi/N)*cot((k-j)*pi/N)
        D = - .5 * (-1).^rif .*D ;              % D(k,j) = -6/12 (-1)^(k-j)*csc((k-j)*pi/N)*cot((k-j)*pi/N)
        D(1:(N+2):end) = - (N^2 + 1 ) /12 ;     % change elements on the diag 
    end
end
