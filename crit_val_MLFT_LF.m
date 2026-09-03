function cv = crit_val_MLFT_LF(L,idx_elect,alpha)
% Estimation of MLFT-LF critical values using Monte Carlo simulations.
% It is considered o conjunto de eletrodos na seguinte ordem:
% FCz, F4, T6, P4, T4, Oz, C4, T5, F7, P3, F3, T3, C3, Fz, Pz, and Cz.
%
% Inputs
%   L: number of nearest-neighbor frequencies
%   N: number of channels
%   alpha: significance level

% Initialization
rng(1234,'twister');
N = length(idx_elect);
Nt = 16;
Nit = 1e6;
F = zeros(Nit,1);

% Neighboring electrodes
Si = {[14,16];
    [7,14];
    [4,5];
    [3,7,15];
    [3,7];
    [15];
    [2,4,5,16];
    [10,12];
    [11,12];
    [8,13,15];
    [9,13,14];
    [8,9,10];
    [10,11,12,16];
    [2,11,16];
    [4,6,10,16];
    [7,13,14,15]};

% Monte Carlo simulation
for I = 1:Nit
    % Generating random values
    Yf0 = randn(1,Nt)+1j*randn(1,Nt);
    YfL = randn(L,Nt)+1j*randn(L,Nt);

    % Applying Laplacian filter
    Yf0_LF = zeros(1,N);
    YfL_LF = zeros(L,N);
    for n = 1:N
        idx = idx_elect(n);
        idx_nn = Si{idx};
        Yf0_LF(:,n) = Yf0(:,idx)-0.25*sum(Yf0(:,idx_nn),2);
        YfL_LF(:,n) = YfL(:,idx)-0.25*sum(YfL(:,idx_nn),2);
    end

    % MLFT estimation
    F(I) = sum(abs(Yf0_LF).^2)/mean(sum(abs(YfL_LF).^2,2));
end

% 95th quantile of F
cv = quantile(F,1-alpha);