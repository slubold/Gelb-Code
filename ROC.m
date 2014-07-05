%% Receiver Operating Characteristic (ROC) curve
% This script plots the ROC for the Gibbs concentration factor 
% and a 5-point detector.

close all; clear all; clc

%% Initialization 
% Define the Fourier modes
nModes = 32;
fourModes = (-nModes:nModes).';

% Generate a (periodic) physical space grid
nGridPts = 256;
xl = -pi; xr = pi; h = (xr-xl)/nGridPts;
x = xl + h*(0:nGridPts-1).';
stencil = h*(-1:1).';

%% Concentration factor definitions
% We use the trigonometric/Gibbs concentration factor
Sipi = 1.85193705198247;					% Value of Si(pi)
cfac_t = pi * sin( pi*abs(fourModes)/nModes )/Sipi;
D = abs(fourModes)/nModes; % This is used to make the computation of cface easier.
p = 1; % p is the order.
cfac_p = p*pi*(abs(fourModes)/nModes).^p; %Polynomial concentration factor. 


%% Noise characteristics
% We consider zero mean, additive white Gaussian nose of variance rho^2
rho2 = 3.5/nModes^2;

%% Covariance Structure
% Compute covariance matrix for the choice of concentration factor and 
% detection stencil
% cov = eye( length(stencil) );
%cov = covariance(stencil, cfac, rho2, nModes);

N = nModes;
k = (-N:N).';						% Fourier modes
len = length(stencil);
cov = zeros(len);
% Compute the covariance matrix.
cov(1,1) = sum((rho2) * real(cfac_t.^2));
cov(2,2) = sum((rho2) * real(cfac_p.^2));
cov(3,3) = sum((rho2) * real(cfac_t.^2));
cov(1,2) = sum((rho2) * real((cfac_t.*cfac_p).*exp(1i.*k.* (stencil(1)-stencil(2)))));
cov(1,3) = sum((rho2) * real((cfac_t.^2).*exp(1i.*k.* (stencil(1)-stencil(3)))));
cov(2,1) = sum((rho2) * real((cfac_t.*cfac_p).*exp(1i.*k.* (stencil(2)-stencil(1)))));
cov(2,3) = sum((rho2) * real((cfac_t.*cfac_p).*exp(1i.*k.* (stencil(2)-stencil(3)))));
cov(3,1) = sum((rho2) * real((cfac_t.^2).*exp(1i.*k.* (stencil(3)-stencil(1)))));
cov(3,2) = sum((rho2) * real((cfac_t.*cfac_p).*exp(1i.*k.* (stencil(3)-stencil(2)))));
% SNR
snr = real( template.'*inv(cov)*template );


%% Function definition
% The template function is a periodic ramp function with unit jump at x = 0
fx = ( (pi-x) .* (x>=0) + (-x-pi) .* (x<0) )/(2*pi);
plot(x,fx);

% Define detector stencil
stencil = h*(-2:2).';

fourMat = exp( 1i*stencil*fourModes.' );

% Now, we generate the Fourier modes of the template response/waveform
fHat = 1./(2i*pi*fourModes); fHat(fourModes==0) = 0;

% Template or jump response
jmpFncCfs = 1i * fHat .* sign(fourModes) .* cfac;
template = real( fourMat*jmpFncCfs );

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

    % Compute the concentration jump approximation
    % For the null hypothesis
    jmpFncCfs = 1i * fHat_null .* sign(fourModes) .* cfac;
    f_null = fourMat * jmpFncCfs;

    % Did the detector register a false alarm?
    pfa = pfa + ( template.'*inv(cov)*f_null > gamma );
	
    % And now, for the alternate hypothesis
    jmpFncCfs = 1i * fHat_alt .* sign(fourModes) .* cfac;
    f_alt = fourMat * jmpFncCfs;

    % Did the detector see the jump?
    pd = pd + ( template.'*inv(cov)*f_alt > gamma );
    	
end

pfa = pfa/nIter; pd = pd/nIter;

%% Plot ROC cruve
figure; plot(pfa, pd, 'r-', 'linewidth', 2); hold on
title 'ROC Curve, exponential concentration factor'
grid on
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