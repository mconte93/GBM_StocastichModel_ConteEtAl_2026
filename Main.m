%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This script simulates the temporal evolution of a glioma growth and 
%   therapy model over a specified number of stochastic realizations. It 
%   calculates finite element updates, applies radiation therapy schedules, 
%   and tracks tumor volume metrics over time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

% --- Configuration and Parameters ---
Parameters; % Load model parameters script
week_of_therapy = 6;
day_of_therapy = 5 * week_of_therapy;
th_T2 = 0.16;  % T2 detection threshold

% --- Load Mesh and Matrix Data ---
load('Mesh_2E_mod.mat');
C = C_mod; 
L = L_mod;
N = size(C, 2);
s = size(L);

load('Matrices_ModStoc.mat', 'M', 'S', 'Pix');
S_M = acc .* S;

load('ABU02_ODF_output_noUnit.mat');
C_new_ODF(1:2, :) = C_new_ODF(1:2, :);
C_new_ODF(3:5, :) = acc .* C_new_ODF(3:5, :);

load('Fiber_2.mat');
Q = Q_2;
Q_mesh = matchQ(Q, C, Pix, h_vox);

% --- Stochastic Simulation Setup ---
SET_of_SIM = 1;
NS = 25;   % Number of stochastic simulations

for ns = (SET_of_SIM * NS - NS) + 1 : SET_of_SIM * NS
    
    InCon_Ter; % Initialize conditions for therapy
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Definition of Time Steps and Horizons
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt = 1;                      % Time step (1 day)
    Tmax = 1 + 7 * week_of_therapy; % Maximum time horizon for weeks
    Tmax1 = 7 * 14;              % Extended simulation horizon
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initial Data and Vector Allocations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = 1;
    step = 1;
    X = U(:, 1);
    
    % --- Therapy Vector Generation (Weekends off) ---
    TTr = zeros(1, Tmax1);
    j = 1;
    for i = 1 : Tmax - 1
        if mod(j, 7) == 0
            TTr(i) = 0;
            TTr(i - 1) = 0;
        else
            TTr(i) = 1;
        end
        j = j + 1;
    end
    TTr = [0, TTr];
    dr = 60 / day_of_therapy;
    
    % --- Stochastic Setup (Ornstein-Uhlenbeck Process) ---
    xi(1) = rand(1, 1);
    mu_R = 1;
    sigma_R = 0.65;
    rng("shuffle");
    for i = 1 : Tmax1
        xi(i + 1) = xi(i) + (mu_R - xi(i)) * dt + sqrt(dt) * sigma_R * xi(i) * randn;
        if xi(i + 1) < 0
            xi(i + 1) = 0;
        end
    end
    
    % --- Metrics Initialization ---
    tau_M_2perc = 0;
    tau_Q = 0;
    timeRec_M = 0;
    
    Tvol(1) = T0;
    Tvol_red(1) = 0;
    
    Tvol_T2(1) = T0_thT2;
    Tvol_T2_red(1) = 0;
    
    Qvol(1) = Q0;
    Qvol_T2(1) = Q0_thT2;
    Qvol_red(1) = 0;
    Qvol_T2_red(1) = 0;
    
    par = [lambda0, k1plus, k2plus, kminus, lambda1, lambda0, KM, speed_M, ...
           mu_M, acc, alfaQ, alfaM, betaQ, betaM, var_eps, r_Q, th_T2, NS, ...
           mu_R, sigma_R];
       
    tic;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main Time-Stepping Loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while T < Tmax1
          
          step = step + 1;
          flag_radio = TTr(step);
          
          % --- Radiation Response Calculation ---
          if TTr(step) == 1
              SM = exp(-alfaM * dr - betaM * dr^2);
              R_M = 1 - SM;
              SQ = exp(-alfaQ * dr - betaQ * dr^2);
              R_Q = 1 - SM;
          else
              CloseZero = 0;
              for k = 1 : step - 1
                  if TTr(k) == 1
                      CloseZero = k;
                  end
              end
              dist = step - CloseZero;
              SM = exp(-alfaM * dr - betaM * dr^2);
              R_M = (1 - SM) * exp(-(dist) / var_eps);
              SQ = exp(-alfaQ * dr - betaQ * dr^2);
              R_Q = (1 - SQ) * exp(-(dist) / var_eps);
          end
          
          L_M = xi(step) * R_M .* ones(N, 1);
          L_Q = xi(step) * R_Q .* ones(N, 1);
        
          P_M = mu_M .* (1 - (U(:, step-1) - Q_evo(:, step-1)) / 2) - L_M;
          
          B_QM = k1plus * SQ * Q_evo(:, step-1) + k2plus * SM * U(:, step-1) + kminus;
          F_Q = lambda1 .* kminus .* k1plus .* SQ ./ (B_QM .^ 2 .* (B_QM + lambda0));
          F_M = lambda1 .* kminus .* k2plus .* SM ./ (B_QM .^ 2 .* (B_QM + lambda0));
          
          [CON_M, b_M] = matrices_M(C, L, C_new_ODF, X, Pix, F_Q, F_M, Q_evo(:, step-1), P_M, h_vox);
       
          % --- Finite Element Method (FEM) Assembly and Solution ---
          A_M = S_M - CON_M;
          A_M = sparse(A_M);
          
          Mbar = M;
          Sbar = A_M;
          Fbar = b_M;
          
          A_FEM = dt * Sbar + Mbar;
          G = Mbar * X + dt * Fbar;
          V = A_FEM \ G;
          
          X = V;
          U(:, step) = max(X, 0);
          
          % --- Normal tissue evolution ---
          dQ_S = -L_Q + r_Q * (1 - (U(:, step) - Q_evo(:, step-1)) / 2) - d_Q .* U(:, step) ./ (1 + U(:, step));
          Q_evo(:, step) = Q_evo(:, step-1) ./ (ones(N, 1) - dt .* dQ_S);
          
          % --- TCP, NTPC, and Control Flag Evaluations ---
          if tau_M_2perc == 0
              tau_M_2perc = tauM_fun(PTV, U(:, step), step, 2);
          end
          
          if tau_Q == 0
              tau_Q = tauQ_fun(L, PTV, Q_evo(:, step), Q_evo(:, 1), step);
          end
          
          % --- Tumor Volume and Radius Computations ---
          [Tvol(step), Tvol_T2(step)] = Tumor_volume_th(C, L, PTV, U(:, step), th_T2);
          
          Tvol_red(step) = ((Tvol(step) - Tvol(1)) ./ Tvol(1));
          Tvol_T2_red(step) = ((Tvol_T2(step) - Tvol_T2(1)) ./ Tvol_T2(1));
          
          % --- Relapse Time Tracking ---
          if T > Tmax
              if timeRec_M == 0 && Tvol(step) >= Tvol(1)
                  timeRec_M = step;
              end
          end
          
          % --- Tissue Volume Metrics ---
          [Qvol(step), Qvol_T2(step)] = Tumor_volume_th(C, L, PTV, Q_evo(:, step), th_T2);
          Qvol_red(step) = ((Qvol(step) - Qvol(1)) ./ Qvol(1));
          Qvol_T2_red(step) = ((Qvol_T2(step) - Qvol_T2(1)) ./ Qvol_T2(1));
          
          T = T + dt;
          
          clear b_M F_Q F_M B_QM P_M; 
    end
    
    toc; 
    
    % --- Save Realization Results ---
    output_filename = ['Results_15days/Evo15day_Num_', num2str(ns), '.mat'];
    save(output_filename, 'PTV', 'U', 'Q_evo', 'par', ...
         'tau_M_2perc', 'tau_Q', 'Tvol', 'Tvol_T2', 'Qvol', 'Qvol_T2', ...
         'Tvol_red', 'Tvol_T2_red', 'Qvol_red', ...
         'Qvol_T2_red', 'timeRec_M', 'NS', 'Tmax1', 'Tmax', ...
         'day_of_therapy', 'week_of_therapy', 'SET_of_SIM');
         
    clear U Q_evo X step tau_M tau_Q xi;
end

return;