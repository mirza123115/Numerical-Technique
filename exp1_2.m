%% ======================
%   BASIC COMMANDS & SPECIAL CHARACTERS
%% ======================
clear all
close all
clc

x = 1:2:9;              % semicolon (;) suppresses output → x = [1 3 5 7 9]

x = 1:2:9; y = 2:10;    % two statements on same line (semicolon used)

x = [1 3 5 ...          % ellipsis (...) continues statement to next line
     7 9];              % x = [1 3 5 7 9]

x = [1:2:9];            % colon (:) creates range → x = [1 3 5 7 9]

x = [1:2:9]';           % apostrophe (') transpose → column vector [1;3;5;7;9]

x = [1:2:9], y = [1:9]  % comma (,) separates commands → two outputs shown

x = 3:3:9;              % x = [3 6 9]
y = x.^2;               % dot (.) operator → element-wise square: [9 36 81]
y = 3*ones(3,1)';       % row vector [3 3 3]
z = x./y;               % element-wise division

%% ======================
%   SIMPLE EXPRESSION
%% ======================
clear all; 
a = 5^2;                % a = 25
b = 3*pi;               % b = 3π
y = exp(a/b);           % exponential
disp(y)                 % display result

%% ======================
%   HELP COMMANDS
%% ======================
help abs                % open MATLAB help for abs()
help cos                % check cos function documentation

%% ======================
%   EXERCISE 1
%% ======================
y = 3*log(sinh(exp(54/6*pi)));  % compute given expression

%% ======================
%   BUILT-IN FUNCTIONS
%% ======================
a = abs(-5)         % absolute value → 5
b = cos(pi/3)       % cosine → 0.5
c = sin(pi/2)       % sine → 1
d = exp(1)          % exponential e^1
e = log(10)         % natural logarithm
f = real(3+4i)      % real part of complex number → 3
g = sqrt(16)        % square root → 4
h = floor(3.9)      % greatest integer ≤ 3.9 → 3
i = ceil(3.1)       % smallest integer ≥ 3.1 → 4

%% ======================
%   MATRICES
%% ======================
A = [1 2 3; 4 5 6; 7 6 3];   % 3x3 matrix
A(1,3)                       % element at row 1 col 3
Z = zeros(2,3);              % 2x3 zero matrix
O = ones(2,3);               % 2x3 ones matrix
I = eye(3);                  % 3x3 identity matrix
A = 0:0.3:3;                 % vector from 0 to 3 with step 0.3
help size; help prod; help length

% EXERCISE 2
A = [1 2 3; 4 5 6; 7 6 54; 65 23 45];
B = 7:1:13.5;
size(A), length(A)           % size and length of A
size(B), length(B)           % size and length of B
A(1:2,2:3)                   % submatrix rows 1-2, cols 2-3
A([1 2],[2 3])               % another submatrix selection
A(1,1) = sin(5);             % assign new value to element

%% ======================
%   MATRIX OPERATIONS
%% ======================
x = 3:3:9;                   % [3 6 9]
y = 3*ones(3,1)';            % [3 3 3]
z = x./y;                    % element-wise division

% EXERCISE 3
A = [2 3; 4 5];
B = [3 4; 6 7];
A+B                          % addition
A*B                          % matrix multiplication
A.*B                         % element-wise multiplication
A/B                          % right matrix division
A\B                          % left matrix division
A.^2                         % element-wise power
A./B                         % element-wise division

%% ======================
%   GRAPHICS
%% ======================
x = 0:0.1:pi; 
y = cos(x);
plot(y);                     % simple plot
plot(x,cos(x),'r');          % red cosine curve
plot(x,y,x,y.^2);            % multiple plots in same figure

% Example with hold
x = 0:0.1:pi; 
y = cos(x);
plot(y); 
hold on; 
plot(x,cos(x),'r');

% Subplot example
x = linspace(0,7); 
y = exp(x);
subplot(2,1,1), plot(x,y);
subplot(2,1,2), semilogy(x,y);

% Bar graph
x = magic(3); 
bar(x); grid on;

% EXERCISE 4: Multiple sin curves
x = 0:0.01:pi;
y1 = sin(x); y2 = sin(2*x); y3 = sin(3*x); y4 = sin(4*x);
plot(x,y1,x,y2,x,y3,x,y4);

%% ======================
%   LOOPS & CONDITIONALS
%% ======================

% Example 6: IF
x = -3;
if x > 0
    a = 10;
elseif x < 0
    a = 11;
elseif x == 0
    a = 12;
else
    a = 14;
end
a

% Example 7: WHILE
x = -10;
while x < 0
    x = x + 1;
end
x

% Example 8: FOR
x = 0;
for i = 1:10
    x = x + 1;
end
x

% Example 9: BREAK
x = -10;
while x < 0
    x = x + 2;
    if x == -2
        break;
    end
end
x

%% ======================
%   USER-DEFINED FUNCTION
%% ======================
% File: cal_pow.m
% function y = cal_pow(x)
%     y = 1 + x^2;
% end

% Main script
% clear all;
% x = 0:1:3;
% t = length(x);
% for i = 1:t
%     val(i) = cal_pow(x(i));
% end
% plot(x,val);

%% ======================
%   EXERCISE 5: VARIANCE
%% ======================
x = 1:1000;
sigma1 = var(x);             % variance using built-in
x = rand(1,1000);            
sigma2 = var(x);             % variance of random data

%% ======================
%   EXERCISE 6: MATRICES & SYSTEM OF EQUATIONS
%% ======================
A = [17 2 3 4; 
     5 6 7 8; 
     9 10 11 12; 
     13 14 15 16];

B = [2 3 4 5; 
     6 7 8 9; 
     10 11 12 13; 
     14 15 16 17];

C = [1 2 3; 
     4 5 6; 
     7 8 9];

y = [4 3 2 1]';              % column vector

AB = A*B; 
BA = B*A;                    % not commutative
% AC = A*C;                  % ERROR: dimension mismatch


% Solve system A*X = b

A = [17 2 3 4;
     5 6 7 8;
     9 10 11 12;
     13 14 15 16];

b = [4; 3; 2; 1];

X = A\b   % solve system A*X = b



% Kaiser Window Based FIR Low-pass Filter Design
% ==============================================
clear; clc; close all;

%% Step 1: Define Inputs
Fs = 1000;          % Sampling frequency (Hz)
Fc = 250;           % Cutoff frequency (Hz)
DeltaF = 50;        % Transition width (Hz)
delta = 0.001;      % Passband ripple
Lmax = 32;          % Series limit (for optional custom Bessel)

%% Step 2: Derived Parameters
fc = Fc/Fs;                     % Normalized cutoff frequency (0 to 0.5)
Dw = 2*pi*DeltaF/Fs;            % Normalized transition width (radians/sample)

A = -20*log10(delta);           % Stopband attenuation (dB)
N = floor((A-8)/(2.285*Dw)) + 1; % Kaiser window length
alpha = (N-1)/2;                 % For linear phase

% Kaiser window beta parameter
if A <= 21
    beta = 0;
elseif A > 21 && A <= 50
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

%% Step 3: Kaiser Window using besseli
n = 0:N-1;
wn = besseli(0, beta*sqrt(1 - ((n-alpha)/alpha).^2)) / besseli(0, beta);

%% Step 4: Ideal Low-pass Filter Impulse Response
hd = zeros(1,N);
for i = 1:N
    if (n(i)-alpha)==0
        hd(i) = 2*fc; % sinc(0) = 1
    else
        hd(i) = sin(2*pi*fc*(n(i)-alpha)) / (pi*(n(i)-alpha));
    end
end

%% Step 5: FIR Coefficients
hFIR = hd .* wn;

%% Step 6: Plot Impulse Response
figure;
stem(n,hFIR,'filled'); grid on;
xlabel('n'); ylabel('h[n]');
title('Impulse Response of FIR Filter (Kaiser Window)');

%% Step 7: Frequency Response using FFT
NFFT = 1024;
H = fft(hFIR, NFFT);
faxis = (0:NFFT-1)/NFFT*Fs;

magdB = 20*log10(abs(H));
phaseDeg = unwrap(angle(H))*180/pi;

figure;
subplot(2,1,1);
plot(faxis(1:NFFT/2+1), magdB(1:NFFT/2+1)); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Magnitude Response');

subplot(2,1,2);
plot(faxis(1:NFFT/2+1), phaseDeg(1:NFFT/2+1)); grid on;
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
title('Phase Response');

%% Step 8: Verification using freqz
figure;
freqz(hFIR, 1, 1024, Fs);
title('Frequency Response using freqz()');
