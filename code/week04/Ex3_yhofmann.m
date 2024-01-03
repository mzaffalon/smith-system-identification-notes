close all
clear all
clc
%% Test 1
test_Signal = randn(1,100);
res_fft = fft(test_Signal,100);
res_self= fourierTransformEx3(test_Signal,100);
if(res_self ~= res_fft)
    disp('False')
else
    disp('Is correct')
end

dif = sum(abs(res_self-res_fft))
%is correct now

%%
% a)
e = randn(1,1024);
mean(e);
var(e);

%%
% b)
res = (abs(fft(e)).^2)/1024;

%[pxx,f] = periodogram(e); Do not use the built in function
%freq = linspace(0,pi,1024); Don't use linspace
freq = 2*pi/1024*(0:1024-1);  %Use this
idx = find(freq>=0 & freq <pi);%Necessary!!!


figure(1)
loglog(freq(idx),res(idx))
xlim([0,pi])


%%
% c)
T = 1;
z = tf('z', T);
P = (z+0.5)/(0.5+ (z+0.5)*(z-0.5)^2);
time_steps = (0:length(e)-1);
w = lsim(P,e,time_steps);

%%
% d)
res_1024 = ((abs(fft(w)).^2)/1024)';

%Questions: Where now only take the ones up to pi?
%Plot goes further to the right than solutions, is solution only axis
%limits different or 

%%
%e)
figure(2)
loglog(freq(idx),res_1024(idx));
hold on

[H, w] = freqresp(P, freq);
magnitude = squeeze(abs(H));
loglog(w(idx),magnitude(idx).^2);
hold on

dif = res_1024' - magnitude;
loglog(freq(idx),abs(dif(idx)));
legend('Periodogram','Plant', 'Difference')
