%% Receiver Operating Characteristic (ROC) curve
% This script plots the ROC for the Gibbs concentration factor 
% and a 5-point detector.

close all; clear all; clc

%% Initialization 
% Define the Fourier modes
nModes = 64;
fourModes = (-nModes:nModes).';

% Generate a (periodic) physical space grid
nGridPts = 256;
xl = -pi; xr = pi; h = (xr-xl)/nGridPts;
x = xl + h*(0:nGridPts-1).';
%% Concentration factor definitions
% Define detector stencil
delta = 1/(2*nModes + 1);
stencil = h*(-1:1).';
fourMat = exp( 1i*stencil*fourModes.' );
D = abs(fourModes)/nModes; % This is used to make the computation of cface easier.
%exponential concentration factor ("cface")
alpha = 5; %order
normal = @(x) (exp(1./(alpha*x.*(x-1))));%the function that goes into the integral of the normalizing constant
C = (pi)./(integral(normal,(2/nModes),1-(2/nModes))); %normalizing constant
cfac_e = C.*D.*exp((1)./(alpha.*D.*(D-1)));
cfac_e(1) = 0; cfac_e(length(cfac_e)) = 0;
%trig concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi), the sine integral                      
cfac_t = pi * sin(pi*abs(fourModes)/nModes)/Sipi; 
%polynomial concentration factor
p = 1; %order of polynomial
cfac_p = p*pi*(abs(fourModes)/nModes).^p; 
figure;
plot(fourModes, cfac_e, fourModes, cfac_t, fourModes, cfac_p);
grid on
title('Concentration factors plotted in Fourier space');
xlabel('k')
ylabel('\sigma(|k|/N)')
legend('Exponential','Trigonometric','Polynomial')
%Here we will call the conc design functions to use other concentration
%factors that are designed by cvx.
cfas = conc_design_Sparse_Jump;  
cfaho = conc_design_Higher_Order; 
% figure 
% plot(fourModes, cfaho); grid on
% title 'Higher-order kernel constraint concentration factor'
cfp4 = conc_design_Problem_4;
% figure 
% plot(fourModes, cfp4); grid on
% title 'problem four concentration factor'


%% Noise characteristics
% We consider zero mean, additive white Gaussian nose of variance rho^2
rho2 = 3.5/nModes^2;

%% Define the covariance matrix here.
N = nModes;
k = (-N:N).';						% Fourier modes
len = length(stencil);
cov = zeros(len);
% Compute the covariance matrix.
cov(1,1) = sum((rho2/2) * real(cfas.^2));
cov(2,2) = sum((rho2/2) * real(cfas.^2));
cov(3,3) = sum((rho2/2) * real(cfas.^2));
cov(1,2) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(1)-stencil(2)))));
cov(1,3) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(1)-stencil(3)))));
cov(2,1) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(2)-stencil(1)))));
cov(2,3) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(2)-stencil(3)))));
cov(3,1) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(3)-stencil(1)))));
cov(3,2) = sum((rho2/2) * real((cfas.^2).*exp(1i.*k.* (stencil(3)-stencil(2)))));
%% Function definition
% The template function is a periodic ramp function with unit jump at x = 0
fx = ( (pi-x) .* (x>=0) + (-x-pi) .* (x<0) )/(2*pi);
plot(x,fx);
grid on
xlabel('x')
ylabel('f(x)')
legend('Sawtooth function')
title('Sawtooth function plot')

% Now, we generate the Fourier modes of the template response/waveform
fHat = 1./(2i*pi*fourModes); fHat(fourModes==0) = 0;

% Template or jump response. Calculate using three different concentration
% factors.
Exponent_factor_1 = exp( 1i*stencil(1) * fourModes.');
Exponent_factor_2 = exp( 1i*stencil(2) * fourModes.');
Exponent_factor_3 = exp( 1i*stencil(3) * fourModes.');
jmpFncCfs_1 = 1i * fHat .* sign(fourModes) .* cfas;
jmpFncCfs_2 = 1i * fHat .* sign(fourModes) .* cfas;
jmpFncCfs_3 = 1i * fHat .* sign(fourModes) .* cfas;
template_1 = real(Exponent_factor_1 * jmpFncCfs_1);
template_2 = real(Exponent_factor_2 * jmpFncCfs_2);
template_3 = real(Exponent_factor_3 * jmpFncCfs_3);
template(1) = template_1; template(2) = template_2; template(3) = template_3; template = template';
% This assigns the values the jump approximation to an element in template, which is the M vector.
snr = template'*inv(cov)*template;
%% Generate ROC curve
% Generate ROC curve using Monte-Carlo methods
nIter = 10000;                      % No. of iterations
nROCPts = 2000;                      % How many data points do we want in 
                                    % ROC curve?
gamma = linspace(-10,10,nROCPts).';	% Range of thresholds to generate ROC

pfa = zeros(nROCPts, 1);
pd  = zeros(nROCPts, 1);

for l = 1:nIter
    
    % Generate realization of noise
    % ZMAWGN
	n = sqrt(rho2/2)*randn(2*nModes+1,1) + ...
                1i*sqrt(rho2/2)*randn(2*nModes+1,1);
    
    % Generate Fourier modes of the null and alternate hypothesis by adding
    % noise
    fHat_null = n;
    fHat_alt  = fHat + n;
    %fHat_alt = template + n;

    % Compute the concentration jump approximation
    % For the null hypothesis
    jmpFncCfs_1 = 1i * fHat_null .* sign(fourModes) .* cfas;
    jmpFncCfs_2 = 1i * fHat_null .* sign(fourModes) .* cfas;
    jmpFncCfs_3 = 1i * fHat_null .* sign(fourModes) .* cfas;
    f_null(1) = Exponent_factor_1 * jmpFncCfs_1;
    f_null(2) = Exponent_factor_2 * jmpFncCfs_2;
    f_null(3) = Exponent_factor_3 * jmpFncCfs_3;

    % Did the detector register a false alarm?
    pfa = pfa + ( template.'*inv(cov)*f_null.' > gamma );
	
    % And now, for the alternate hypothesis
    jmpFncCfs_1 = 1i * fHat_alt .* sign(fourModes) .* cfas;
    jmpFncCfs_2 = 1i * fHat_alt .* sign(fourModes) .* cfas;
    jmpFncCfs_3 = 1i * fHat_alt .* sign(fourModes) .* cfas;
    f_alt(1) = Exponent_factor_1 * jmpFncCfs_1;
    f_alt(2) = Exponent_factor_2 * jmpFncCfs_2;
    f_alt(3) = Exponent_factor_3 * jmpFncCfs_3;

    % Did the detector see the jump?
    pd = pd + ( template.'*inv(cov)*f_alt.' > gamma );
    	
end
pfa = pfa/nIter; pd = pd/nIter;

%% Plot ROC cruve
figure; plot(pfa, pd, 'r-', 'linewidth', 2); hold on
grid on
title 'ROC Curve'
xlabel 'Probability of False Alarm, Pfa'; 
ylabel 'Probability of Detection, Pd'

% Analytic ROC curve 
% This can be worked out in terms of the Q (right tail Gaussian 
% probability) function; see, for example, Fundamentals of Statistical 
% Signal Processing: Detection Theory (Vol 1), Steven Kay
pfa = linspace(0, 1, nROCPts);
pd = Q(Qinv(pfa) - sqrt(snr));
plot(pfa, pd, 'b--', 'linewidth', 2)
legend('Monte Carlo Simulation', 'Analytic')