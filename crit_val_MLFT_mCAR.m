function CV = crit_val_MLFT_mCAR(L,N,alpha)
% Estimation of MLFT-mCAR critical values using Monte Carlo simulations.
% 
% Inputs
% L: number of nearest-neighbor frequencies
% N: number of channels
% alpha: significance level


% Initialization
rng(1234,'twister'); 
Nit = 1e6;
F = zeros(Nit,1);

% Monte Carlo simulation
parfor I = 1:Nit
    % Generating random values
    Yf0 = randn(1,N)+1j*randn(1,N);
    YfL = randn(L,N)+1j*randn(L,N);

    % Applying mCAR
    Yf0_mCAR = Yf0-mean(Yf0,2);
    YfL_mCAR = YfL-mean(YfL,2);

    % MLFT estimation
    F(I) = sum(abs(Yf0_mCAR).^2)/mean(sum(abs(YfL_mCAR).^2,2));
end

% 95th quantile of F
CV = quantile(F,1-alpha);