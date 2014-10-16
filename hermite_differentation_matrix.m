function D2 = hermite_differentation_matrix(N)
    
    
    % calculate the hermite polynomials zero
    T = 1/sqrt(2) * eig(diag(sqrt((1:N-1)),1) + diag(sqrt((1:N-1)),-1));
    
    % utils: calculate T(i) - T(j) for all i,j
    minus = T*ones(1,N) - ones(N,1) * (T');
    minusplus = minus+eye(N);
    % pi1T(j) = pi'(t_j)
    pi1T = (prod((minusplus)'))';
    
    % pi2T(j) = pi''(t_j)
    pi2T = zeros(N,1);
    for i = (1:N)
        mask = ones(N,1);
        mask(i) = 0;
        minusplusMOD = minusplus;
        minusplusMOD(:,i)= mask;
        pi2T = pi2T + 2*(prod((minusplusMOD)'))';
    end
    
    % pi3T(j) = pi'''(t_j)
    pi3T = zeros(N,1);
    for i = (1:N)
        for j = (1:N)
            if i ~= j
                maski = ones(N,1);
                maski(i) = 0;
                maski(j) = 0;
                maskj = ones(N,1);
                maskj(i) = 0;
                maskj(j) = 0;
                minusplusMOD = minusplus;
                minusplusMOD(:,i)= maski;
                minusplusMOD(:,j)= maskj;
                pi3T = pi3T + 3*(prod((minusplusMOD)'))';
            end
        end
    end
    
    
    % calculate the first order derivative matrix
    D1 = (pi1T * ones(1,N)) ./ (ones(N,1) * pi1T') ./minusplus ;
    D1 = D1 .* (ones(N) - eye(N));
    D1 = D1 + diag(0.5 * pi2T ./ pi1T );
    
    % calculate the second order derivative matrix
    D2 =  ((pi2T * ones(1,N)) ./ (ones(N,1) * pi1T') - 2*D1 ) ./minusplus;
    D2 = D2 .* (ones(N) - eye(N));
    D2 = D2 + diag(1/3 * pi3T ./ pi1T);


