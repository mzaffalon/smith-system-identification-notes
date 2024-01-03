clear all
close all
%%1) 
% Define the model parameters
a = 5/8;
b = 11/10;
c = 5/6;
covariance = 0.05;
G = tf([b, c], [1, -a, 0], 1);

% Define the number of periods and time steps
num_periods = 2;
tau_max = 80;
num_samples = 100;  %Magic :)
num_sinusiods = tau_max/2; 

% Initialize the input signal u(k)
u = generateSinusoidalInput(num_samples, num_periods, num_sinusiods);


% No noise
    y_no_noise = generateOutput(u,zeros(num_samples*num_periods,1),a,b,c,num_periods,num_samples,G);

    g0 = timeEstimation(u,tau_max,y_no_noise);

    g = zeros(tau_max+1,3);

% Simulate the system and estimate the impulse response
for realization = 1:3
    % Generate random noise
    v = sqrt(covariance) * randn(num_samples*num_periods,1);
   
    y = generateOutput(u,v,a,b,c,num_periods,num_samples,G);
    
    g(:,realization) = timeEstimation(u,tau_max,y);
end

% Calculate Average 
g_hat_final = mean(g,2);

error1 = vecnorm(g_hat_final - g0,2,1);



%%
% 2)
max = 80;
g = zeros(tau_max+1,max-1);
for cur_period = 2:max
    u = generateSinusoidalInput(num_samples, cur_period, num_sinusiods);

    v = sqrt(covariance) * randn(num_samples*cur_period,1);
   
    y = generateOutput(u,v,a,b,c,cur_period,num_samples,G);
    
    g(:,cur_period-1) = timeEstimation(u,tau_max,y);  
end

error = vecnorm(g - g0,2,1);
figure(1)
loglog(2:max,error);
title('Error for part 2')

%%
% 3)
g_rand = zeros(tau_max+1,max-1);
for cur_period = 2:max
    u = 4*(rand(cur_period*num_samples,1) - 0.5);

    v = sqrt(covariance) * randn(num_samples*cur_period,1);

    y = generateOutput(u,v,a,b,c,cur_period,num_samples,G);
    
    g_rand(:,cur_period-1) = timeEstimation(u,tau_max,y);  
end




error_rand = vecnorm(g_rand - g0,2,1);
figure(2)
plot(2:max,error_rand);
hold on
plot(2:max,error)
legend('Random', 'Sinusoidal')
title('Error for part 3')




%%
% 4)
u = generateSinusoidalInput(num_samples, num_periods, num_sinusiods);
g4 = zeros(tau_max+1,200);
for realization = 1:200
    % Generate random noise
    v = sqrt(covariance) * randn(num_samples*num_periods,1);
   
    y = generateOutput(u,v,a,b,c,num_periods,num_samples,G);
    
    if realization > 1
        g4(:,realization) = (timeEstimation(u,tau_max,y) + g4(:,realization-1));  
    else
        g4(:,realization) = timeEstimation(u,tau_max,y);
    end
end

for i = 1:size(g4,2)
    g4(:,i) = g4(:,i)/i;
end

error4 = vecnorm(g4 - g0,2,1);
figure(3)
loglog(1:200,error4);
hold on
x = 1:200;
loglog(x,1./x)
hold on
loglog(x,1./sqrt(x))
legend('Error of g', '1/N', '1/sqrt(N)')
title('Error for part 4')


%%
%Functions

% Function to generate sinusoidal input signal
function u = generateSinusoidalInput(samples, periods, sinusoids)
 omega = 2*pi/samples*round(samples/2*(1:sinusoids)'/(sinusoids+1)); % new

 thetas = pi*(1:sinusoids)'.*((1:sinusoids)'-1)/sinusoids; % Schroeder phase


 k = (0:samples-1)';
 u = zeros(samples,1);
 for i = 1:sinusoids
     u = u + 2*cos(k*omega(i)+thetas(i));
 end

 u = 2*u /norm(u,"inf");

 u = repmat(u,periods,1);


end


function y = generateOutput(u,noise,a,b,c, periods, samples, G)
    y = zeros(periods*samples,1);
    y_neg1 = 0;
    u_neg1 = 0;
    u_neg2 = 0;
    y(1,1) = a*y_neg1 + b*u_neg1 + c*u_neg2 + noise(1,1);
    y(2,1) = a*y(1) +b*u(1) +c*u_neg2 + noise(2,1);
   y_neg1 = y(2);
   u_neg1 = u(2);
   u_neg2 = u(1);
    for i = 3:length(u)
        y(i) = a*y_neg1 +b*u_neg1 +c*u_neg2;
        y_neg1 = y(i);
        u_neg1 = u(i);
        u_neg2 = u(i-1);
    end

    % y = lsim(G,u); % Better this than the above one

    y = y + noise;



end


function g = timeEstimation(u, tau_max, y)
    c_toep = u;
    r = [u(1), zeros(1,tau_max)];
    phi = toeplitz(c_toep,r);
    
    g = phi\y;
    % other approach
    % phy2 = zeros(length(u), tau_max);
    % for i = 1:tau_max+1
    %     phy2(i:end,i) = u(1:end-i+1);
    % end
    % g2 = phy2\y;

end
