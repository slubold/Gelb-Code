%%---------------------------------------------------
% Iterative routine to design concentration factors
% needs CVX to be installed
%%---------------------------------------------------

function sig = conc_design_Problem_4

%% Parameters
% No. of coefficients
N = 64;
% these are the Fourier modes
k = -N:N;

% Reconstruction grid
npts = 257;
x = linspace(-pi,pi,npts).';

% Fourier exponential (DFT) matrix
kern = exp(i*x*k);

% set up some matrices which recreate the form of W_0, W_1, W_2,...
wts = i*sign(k)./(i*k); wts(N+1) = 0; wts = diag(wts);
wts_der = i*sign(k)./((i*k).^2); wts_der(N+1) = 0; wts_der = diag(wts_der);
wts_2der = i*sign(k)./((i*k).^3); wts_2der(N+1) = 0; wts_2der = diag(wts_2der);

% for missing data, enter band of missing modes here
% ind = (k>=20).*(k<=30); ind =(ind==1);

%% Iterative design
% this is the bit where we call CVX to solve the variational formulation

cvx_begin
    % conc. factor is an even function
    variable cfac(N);
    
    % need to only solve for k = 1:N
    sig = [flipud(cfac); 0; cfac];
    
    % these are the different kernels W_0, W_1, ....
    % Note: S_N[f] = [f]*W_0 + [f']*W_1 + [f'']*W_2 ....
    
    ker = kern*wts*sig/(2*pi);
    der_ker = kern*wts_der*sig/(2*pi);
    der_2ker = kern*wts_2der*sig/(2*pi);
    
    % choose objective function
    %minimize norm(ker,2)
    minimize norm(ker,2)
    
    % constraints
    subject to
         % normalization for correct jump height
         real(ker(x==0))==1;
         %other possible constraints
         % sig>=0;
         %sig(1)==0;
         norm(ker(abs(x)>=.2),inf)<=1e-4;
     
cvx_end

% 
% %% Plot results
% plot conc. factor
figure; plot(k, sig, 'k-.', 'Linewidth',1); xlabel k; ylabel \sigma(|k|/N)

% 
% % plot response to a test function
% % Test function
% % Ramp
% %fhat = -i./(2*pi*k); fhat(N+1)=0; fhat = fhat.';
% %fx = (x+pi)/2.*(x<=0) + (x-pi)/2.*(x>0); fx = -fx/pi;
% 
% % test function with variation
% fx = sin(6*x).*(x<0) + exp(-3*x).*(x>=0);
% fhat = exp(-i*k.'*x.')*fx(:)/(npts);
% 
% % another test function
% %fx = sin(3*x).*(x<-pi/3) + tanh(x).*(x>=-pi/3).*(x<pi/2) - (3*x/pi-3).*(x>=pi/2);
% %fhat = exp(-i*k.'*x.')*fx(:)/(npts);
% 
% % if you have missing data
% %ind = (k>=20).*(k<=30); ind =(ind==1);
% %fhat(ind) = 0;
% 
% % compute the jump response
% resp = real(kern*(i*sign(k.').*sig.*fhat));
% figure; plot(x, fx, 'k--', 'Linewidth',1); hold on; plot(x, resp, 'k', 'Linewidth',1)
% xlabel x; ylabel [f](x); legend('f','S_N^\sigma[f]'); xlim([-pi pi])