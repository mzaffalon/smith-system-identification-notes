
clear all
close all
clc



%% Ex 3
close all

data1 = importdata("experiment1.dat").data;
t1 = data1(:,1),
y1 =  data1(:,2)
std1= data1(:,3)
reg = [ones(length(y1),1), t1, t1.^2 ];
theta = reg\y1


R = diag(std1.^2);
Z_t_direct = theta /y1; % --> Not good don't do this as probably some numerical stability issues
% or with formula
Z_t_formula = inv(reg'*inv(R)*reg)*reg'*inv(R);



% dont do cov_theta1 = Z_t_direct*R*transpose(Z_t_direct)

% Both work well
cov_theta2 = Z_t_formula*R*Z_t_formula'
%or direct formula
cov_theta3 =inv(reg'*inv(R)*reg);

cov1 = sqrt(cov_theta3(1,1))
cov2 = sqrt(cov_theta3(2,2))
cov3 = 2*sqrt(cov_theta3(3,3))

figure;
errorbar(t1,y1,std1, 'LineStyle','none','Marker','.')
hold on
t_cont = t1(1):0.01:t1(end);

y_pred= reg*theta;
y0 = theta(1,1)
v0 = theta(2,1)
a = theta(3,1)
plot(t_cont, polyval([a v0 y0],t_cont))

g_exp1 = 2*a

hold off

%%
%%Ex 3 part 2
data2 = importdata("experiment2.dat").data;
t2 = data2(:,1),
y2 =  data2(:,2)
std2= data2(:,3)
reg2 = [ones(length(y2),1), t2, t2.^2 ];
%w = cov with formula from ex1
R2 = diag(std2.^2);
cov =inv(reg2'*inv(R2)*reg2);
cov12 = sqrt(cov(1,1))
cov22 = sqrt(cov(2,2))
cov32 = 2*sqrt(cov(3,3))

w = R2;


%or for weighted:
theta2 = inv(reg2'*inv(w)*reg2)*reg2'*inv(w)*y2
%



figure;
errorbar(t2,y2,std2, 'LineStyle','none','Marker','.')
hold on



y02 = theta2(1,1)
v02 = theta2(2,1)
a2 = theta2(3,1)
t_cont = t2(1):0.01:t2(end);
plot(t_cont, polyval([a2 v02 y02],t_cont))

g_exp2 = 2*a2

hold off

%% Ex 4
w0 = 10;
theta0 = pi/8;
data3 = importdata("experiment3.dat").data;
t3 = data3(:,1);
y3 =  data3(:,2);

figure;
plot(t3,y3)

std3= data2(:,3);
reg3 = [ones(size(t3,1),1), t3, 0.5* t3.^2, cos(theta0 + w0*t3)];
theta3 = reg3\y3;

disp(theta3')

true_val = [12.42, 44.19, -6.42, 1.53];
disp("True values")
disp(true_val)

