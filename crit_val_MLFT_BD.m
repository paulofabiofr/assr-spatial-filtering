function CV = crit_val_MLFT_BD(L,N,alpha)
% Estimation of MLFT-BD critical values using Monte Carlo simulations.
% 
% Inputs
%   L: number of nearest-neighbor frequencies
%   N: number of channels
%   alpha: significance level


% Initialization
rng(1234,'twister'); 
Nit = 1e6;
bd = nchoosek(1:N,2); % bipolar devivations
F = zeros(Nit,1);

% Monte Carlo simulation
parfor I = 1:Nit
    % Generating random values
    Yf0 = randn(1,N)+1j*randn(1,N);
    YfL = randn(L,N)+1j*randn(L,N);

    % Applying mCAR
    bdYf0 = Yf0(:,bd(:,1))-Yf0(:,bd(:,2));
    bdYfL = YfL(:,bd(:,1))-YfL(:,bd(:,2));
    Yf0_BD = [Yf0, bdYf0];
    YfL_BD = [YfL, bdYfL];

    % MLFT estimation
    F(I) = sum(abs(Yf0_BD).^2)/mean(sum(abs(YfL_BD).^2,2));
end

% 95th quantile of F
CV = quantile(F,1-alpha);