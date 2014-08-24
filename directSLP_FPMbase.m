function Phi = directSLP_FPMbase(N)
    Xi = 2*pi/N*(0:N) ;
    if mod(N,2)
        Phi = @(xi) 1/N * sin(N/2*(xi - Xi)).*cot((xi - Xi)/2) ; 
    else
        Phi = @(xi) 1/N * sin(N/2*(xi - Xi)).*csc((xi - Xi)/2) ; 
    end
end
