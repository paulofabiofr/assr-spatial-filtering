function cv = crit_val_MLFT_CAR(L,N,Nt,alpha)
% Estimation of MLFT-CAR critical values using Monte Carlo simulations.
% 
% Inputs
%   L: number of nearest-neighbor frequencies
%   N: number of channels
%   Nt: total number of electrodes
%   alpha: significance level

% Initialization
rng(1234,'twister'); 
Nit = 1e6;
F = zeros(Nit,1);

% Monte Carlo simulation
parfor I = 1:Nit
    % Generating random values
    Yf0 = randn(1,Nt)+1j*randn(1,Nt);
    YfL = randn(L,Nt)+1j*randn(L,Nt);

    % Applying CAR
    Yf0_CAR = Yf0(:,1:N)-mean(Yf0,2);
    YfL_CAR = YfL(:,1:N)-mean(YfL,2);

    % MLFT estimation
    F(I) = sum(abs(Yf0_CAR).^2)/mean(sum(abs(YfL_CAR).^2,2));
end

% 95th quantile of F
cv = quantile(F,1-alpha);