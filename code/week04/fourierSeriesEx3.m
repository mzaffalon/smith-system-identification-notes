% function f =fourierSeriesEx3(u, N)
% 
% f = zeros(N,1);
% for n = 1:N
%     w_n = 2*pi *(n-1)/N;
%     temp = 0;
%     for k = 1:n
%         f(n) = temp + u(k)*exp(-1j*w_n*(k-1));
%         temp = f(n);
%     end
% 
% end

% end

function [X] = fourierSeriesEx3(x,N)            
    X = zeros(1, N);            % Initialize the DFT result
    for k = 1:N
        X(k) = 0;
        for n = 1:N
            X(k) = X(k) + x(n) * exp(-1j*2*pi*(k-1)*(n-1)/N);
        end
    end
end
