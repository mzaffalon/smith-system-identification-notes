v = randn(100,1)*sqrt(0.3) + 0.5;
f = @(theta) -sum((y - (theta(1)*x + theta(2)*(2*x.^2 -1) + v)).^2)

theta_initial = [0,0]
options = optimset('Display', 'iter', 'TolFun', 1e-6);
estimated_theta = fminunc(f, theta_initial, options);

fprintf('Estimated theta1: %.4f\n', estimated_theta(1));
fprintf('Estimated theta2: %.4f\n', estimated_theta(2));
%%

% Number 1

for i = 1:size(x,1)
w_k = [x(i),  2*x(i)^(2)-1];
w(i,1) = w_k(1);
w(i,2) = w_k(2);
end


%better 
w = [x 2*x.^2-1];




theta = inv(transpose(w)*w)*transpose(w)*(y-0.5)


% Number 2
var_v = 0.3;
mu_v = 0.5;
var_theta =0.02;
mu_theta = [1.3; 0.9];
theta_MAP = inv((1/var_v)*transpose(w)*w + (1/var_theta)*eye(2)) * ((1/var_v)*transpose(w)*(y-mu_v)+mu_theta/var_theta)

%Number 3
error_ML = 0;
error_MAP =0;
for i = 1:size(x_v,1)
error_ML = error_ML + (theta(1)*x_v(i) + theta(2)*(2*(x_v(i))^2-1) + 0.5 - y_v(i))^2;
error_MAP = error_MAP + (theta_MAP(1)*x_v(i) + theta_MAP(2)*(2*(x_v(i))^2-1) + 0.5 - y_v(i))^2;

end
error_ML =error_ML/size(x_v,1)
error_MAP = error_MAP/size(x_v,1)