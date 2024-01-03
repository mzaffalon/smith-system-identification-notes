function [u, Ru, y, Ruy, Ry,  etfe_raw, spec_raw, etfe_smooth, spec_smooth] = HS2023_SysID_Exercise_05_20914032()

%% Instructions
% Student submissions should be contained in this base script.
% Change the filename and function-name `12345678' to your personal Legi
% number before execution and submission.
%
%% Dependencies
% The following TA-provided files are included in the assignment will be 
% required in order to form solutions:
%
%   HS2023_SysID_Exercise_05_GenerateData.p
%
%% Outputs:
%   u:              The input (chirp) signal applied to the system u(k)
%   Ru:             The autocorrelation of the input (input_u)
%   y:              Output y(k) of the system
%   Ruy:            Cross correlation between u(k) and y(k)
%   Ry:             Autocorrelation of y(k)
%   etfe_raw:       Empirical Transfer Function Estimate
%   spec_raw:       Spectral estimate
%   etfe_smooth:    Smoothed ETFE by Hann window and inverse variance weighting
%   spec_smooth:    Smoothed spectral estimate by Hann window and inverse variance weighting
%% Policies
%
% The only allowable external toolbox to MATLAB is the Control Systems
% Toolbox. MATLAB default commands such as 'fft' or 'fftshift' are allowed.
%
% It is impermissible to other use MATLAB toolboxes, such as the Computer Vision
% Toolbox, the Econometrics Toolbox, the Deep Learning Toolbox, etc.
%
% All student-defined functions should be written at the bottom of this   
% file. The homework code-submission will be this single file.
%
% The code of academic conduct remains in effect when preparing this
% homework.
%
%% Output format specification:
% The output will be stored in the class 'out'. Refer to output_ex_5.m for
% more detail about desired fields and structures


%% Begin routines
% Extract Legi from Filename
name=mfilename;
LegiNumber= name(end-7:end);



% In midterm for plots specify:!!!!!!!!!!!!!!!!!
% Title
% Axis
% Legend

%% problem parameters

%time horizon of the data
T = 1024;

%% responses to problem areas



%% Part 1: Input design

%create functions that will generate inputs under the specified
%requirements      
f0 = 0;           
f1 = 0.3;         
k = (0:T-1);  

%chirp input
u = cos(2*pi*(f0*k+(f1-f0)*k.^2)/(2*T))';
figure(10)

%% Part 2: Input and Autocorrelation


figure(1)
[Ru, valueLag] = intcor(u,u);
stem(valueLag,Ru) % better also use plot
title('Autocorrelation of u')
% from solutions
Ru2 = (ifft(fft(u).*conj(fft(u))))/T; % To shift it to 0 use fftshift
figure(100)
plot(Ru2)

%% Part 3: Gather output

%execute the black-box system with computed inputs to generate noisy
%system responses
y = HS2023_SysID_Exercise_05_GenerateData(20914032, u);



%% Part 4: Output Correlations

%cross correlations between input and output
[Ruy, ~] = intcor(u,y);
% Ryu = conj(Ruy);
% from solution
Ruy2 = ifft(fft(u).*conj(fft(y)))/T;

Ryu = intcor(y, u);


figure(11)
plot(Ruy)
figure(12)
plot(Ruy2)

%autocorrelations of input
[Ry, ~] = intcor(y,y);
Ry2 =  ifft(fft(y).*conj(fft(y)))/T;



%% Part 5: Transfer Function Estimates

%compute empirical transfer function estimates
etfe_raw = fft(y)./fft(u);


spec_raw = fft(Ryu)./fft(Ru); % WATCH OUT is yu and not uy

omega = 2*pi/T *(0:T/2-1);
figure(2)
loglog(omega,abs(spec_raw(1:length(omega),:)));
title('Magnitude of estimated TF from spectral function')
ylim([0.1  10])
xlim([6*10^-3 3])

figure(3)
semilogx(omega,angle(spec_raw(1:length(omega),:)));
title('Phase of estimated TF from spectral function')
ylim([-4 4])


figure(4)
loglog(omega,abs(etfe_raw(1:length(omega),:)));
title('Magnitude of estimated TF from etfe')
ylim([0.1  10])
xlim([6*10^-3 3])

figure(5)
semilogx(omega,angle(etfe_raw(1:length(omega),:)));
title('Phase of estimated TF from etfe')
ylim([-4 4])


%% Part 6: Smoothed Transfer Functions

%smooth the ETFE using windowing and inverse variance weighting
U = fft(u);
[omega6, Wg] = WfHann(20,T); % better gamma = 100
zidx = find(omega6 ==0);
omega6 = [omega6(zidx:T); omega6(1:zidx-1)];
Wg = [transpose(Wg(zidx:T));transpose(Wg(1:zidx-1))];
a = U.*conj(U);
etfe_smooth = 0 *etfe_raw;

for wn = 1:T
    Wnorm = 0;
    for xi = 1:T
        widx = mod(xi-wn,T)+1;
        etfe_smooth(wn) = etfe_smooth(wn)+Wg(widx)*etfe_raw(xi)*a(xi);
        Wnorm = Wnorm +Wg(widx)*a(xi);
    end
    etfe_smooth(wn) = etfe_smooth(wn)/Wnorm;
end

idx = find(omega6>0 & omega6<pi);

%Plots
figure(6)
loglog(omega6(idx), abs(etfe_smooth(idx)))
title('Magnitude of estimated TF from smoothed ETFE')
figure(7)
semilogx(omega6(idx),angle(etfe_smooth(idx,:)));
title('Phase of estimated TF from smoothed ETFE')



%Same for spec
spec_smooth = 0*spec_raw;
for wn = 1:T
    Wnorm = 0;
    for xi = 1:T
        widx = mod(xi-wn,T)+1;
        spec_smooth(wn) = spec_smooth(wn)+Wg(widx)*spec_raw(xi)*a(xi);
        Wnorm = Wnorm +Wg(widx)*a(xi);
    end
    spec_smooth(wn) = spec_smooth(wn)/Wnorm;
end


figure(8)
loglog(omega6(idx), abs(spec_smooth(idx)))
title('Magnitude of estimated TF from smoothed spectral')
figure(9)
semilogx(omega6(idx),angle(spec_smooth(idx,:)));
title('Phase of estimated TF from smoothed spectral')




%% All done!

end

%% Student-defined functions
%[functions go here]
function [R, h] = intcor(u, y)
    % Calculates the correlation between two vectors u, y (periodic signal)

    M = length(u);
    h = -M/2+1:1:M/2;
    R = zeros(M, 1);
    
    for i = 1:M
        for k = 1:M
            if k - h(i) <= 0
                R(i) = R(i) + u(k) * y(k - h(i) + M); % Periodic extension
            else
                if(k - h(i)>M)
                    R(i) = R(i) + u(k) * y(k - h(i)-M);
                
                else
                    R(i) = R(i) + u(k) * y(k - h(i));
                end
            end
        end
    end
    R = R/M;
end

function [omega,WHann] = WfHann(gamma,N)
%------------------------------------------------------------------
%
%   [omega,WHann] = WfHann(gamma,N)
%
%   Create a frequency domain Hann window with width parameter gamma
%   and data length N.  The Hann window is a raised cosine.
%
%   Roy Smith,  18 October, 2017.
%
%                6 November, 2017.  Fixed bug in N even indexing.
%
%------------------------------------------------------------------

if nargin == 0,
    disp('Syntax: [omega,W] = WfHann(gamma,N)')
    return
elseif nargin ~= 2,
    error('incorrect number of input arguments (2 expected)')
    return
end

%   basic parameter checking
if length(gamma) > 1,
    error('Width parameter, gamma, must be a scalar');
end
if round(gamma) ~= gamma,
    error('Width parameter, gamma, must be an integer');
end
if gamma < 1,
    error('Width parameter, gamma, must be positive');
end
if length(N) > 1,
    error('Calculation length, N, must be a scalar');
end
if round(N) ~= N,
    error('Calculation length, N, must be an integer');
end
if N < 1,
    error('Calculation length, N, must be positive');
end

%   The simplest approach is to define the window in the time domain and
%   then transform it to the frequency domain.

lags = [floor(-N/2+1):floor(N/2)]';
wHann = 0*lags;
idx = find(abs(lags) <= gamma);
wHann(idx) = 0.5*(1+cos(pi*lags(idx)/gamma));

%   
zidx = find(lags==0);    % index of the zero point.

wH_raw = fft([wHann(zidx:N);wHann(1:zidx-1)]);
WHann(zidx:N) = wH_raw(1:N-zidx+1);  % shift +ve freq to end
WHann(1:zidx-1) = wH_raw(N-zidx+2:N);% shift > pi freq to beginning
WHann = real(WHann);   % should actually be real
omega = 2*pi/N*lags;

return
  
end
%------------------------------------------------------------------

function [lags,wHann] = WtHann(gamma,N)
%------------------------------------------------------------------
%
%   [lags,wHann] = WtHann(gamma,N)
%
%   Create a Hann window with width parameter gamma and data length N.
%   The Hann window is a raised cosine.
%
%   Roy Smith,  18 October, 2017.
%
%------------------------------------------------------------------

if nargin == 0,
    disp('Syntax: [lags,w] = WtHann(gamma,N)')
    return
elseif nargin ~= 2,
    error('incorrect number of input arguments (2 expected)')
    return
end

%   basic parameter checking
if length(gamma) > 1,
    error('Width parameter, gamma, must be a scalar');
end
if round(gamma) ~= gamma,
    error('Width parameter, gamma, must be an integer');
end
if gamma < 1,
    error('Width parameter, gamma, must be positive');
end
if length(N) > 1,
    error('Calculation length, N, must be a scalar');
end
if round(N) ~= N,
    error('Calculation length, N, must be an integer');
end
if N < 1,
    error('Calculation length, N, must be positive');
end

lags = [floor(-N/2+1):floor(N/2)]';
wHann = 0*lags;
idx = find(abs(lags) <= gamma);
wHann(idx) = 0.5*(1+cos(pi*lags(idx)/gamma));

return
  
end
%------------------------------------------------------------------


