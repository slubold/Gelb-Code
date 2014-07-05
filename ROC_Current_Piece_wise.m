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

%% Concentration factor definitions
% We use the trigonometric/Gibbs concentration factor
% Sipi = 1.85193705198247;					% Value of Si(pi), the sine 
                                            % integral
%cfac = pi * sin(pi*abs(fourModes)/nModes)/Sipi;   

% %polynomial concentration factor
% p = 1; %order of polynomial
% cfac_p = p*pi*(abs(fourModes)/nModes).^p; 
% cfac = cfac_p;

D = abs(fourModes)/nModes; % This is used to make the computation of cface easier.
%exponential concentration factor ("cface")
alpha = 6; %order
normal = @(x) (exp(1./(alpha*x.*(x-1))));%the function that goes into the integral of the normalizing constant
C = (pi)./(integral(normal,(1/nModes),1-(1/nModes))); %normalizing constant
cfac_e = C.*D.*exp((1)./(alpha.*D.*(D-1)));
cfac_e(1) = 0; cfac_e(length(cfac_e)) = 0;
cfac = cfac_e;

%% Function definition
% The template function is a periodic ramp function with unit jump at x = 0
% fx = ( (pi-x) .* (x>=0) + (-x-pi) .* (x<0) )/(2*pi);
fx = (cos(x)).*(x>=0) +  (-pi-x).*(x<0) ; plot(x,fx) 
grid on
title ('test function, [-\pi,\pi]')

%%
% Define detector stencil
%stencil = h*(-1:1).';
%fourMat = exp( 1i*stencil*fourModes.' );

% Now, we generate the Fourier modes of the template response/waveform
fHat = 1./(2i*pi*fourModes); fHat(fourModes==0) = 0;

% Template or jump response
jmpFncCfs = 1i * fHat .* sign(fourModes) .* cfac;
template = real( fourMat*jmpFncCfs );

%% Noise characteristics
% We consider zero mean, additive white Gaussian nose of variance rho^2
rho2 = 3.5/nModes^2;

%% Covariance Structure
% Compute covariance matrix for the choice of concentration factor and 
% detection stencil
% cov = eye( length(stencil) );
cov = covariance(stencil, cfac, rho2, nModes);

% SNR
snr = real( template.'*inv(cov)*template );

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
title 'ROC Curve, trigonometric concentration factor'
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