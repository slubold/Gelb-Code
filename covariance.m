function covr = covariance(stencil, cfac, rho2, N)
%% covariance.m
% Usage: covr = covariance(stencil, cfac, rho2, N) 
% This script computes the covariance matrix associated with the
% concentration method on equispaced modes
% Inputs:	stencil - physical space edge detection stencil
%           cfac    - concentration factor
%           rho2    - noise variance
%           N       - no of Fourier modes used
%
% Outputs:	covr    - covariance matrix

k = (-N:N).';						% Fourier modes
len = length(stencil);
covr = zeros(len);                  % Compute covariance matrix
for a = 1:len
	for b = 1:len
		covr(a,b) = (rho2/2) * real( (cfac.^2).' * ...
                                    exp(1i*k*( stencil(a)-stencil(b) )) ) ;
	end
end

