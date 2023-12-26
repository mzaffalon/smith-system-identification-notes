function [X] = fourierTransformEx3(x,N)            
    X = zeros(1, N);            % Initialize the DFT result
    for k = 1:N
        X(k) = 0;
        for n = 1:N
            X(k) = X(k) + x(n) * exp(-1j*2*pi*(k-1)*(n-1)/N);
        end
    end
end
