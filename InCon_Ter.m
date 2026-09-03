%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Loads pre-computed baseline data from a non-therapy evolution run, 
%   and initializes the numerical tumor and normal tissue populations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Initialization ---
U_0 = zeros(N, 1);    % Preallocate numerical solution for the tumor equation

% --- Load Baseline Data ---
% Load final time-step data from the non-therapy simulation
load('Evo_noTer_th16_05margin.mat', 'U', 'PTV');
U_0 = U(:, end);
PTV_0 = PTV(:, end);

clear U PTV;          

% --- Process Tumor Solution (U) ---
U = U_0;
% Ensure non-negative values for the tumor cell population
for i = 1 : size(U, 1)
   if U(i) < 0
      U(i) = 0;
   end
end

% --- Process PTV and Normal Tissue Population ---
PTV = PTV_0;        
Q_evo = Q_mesh(3, :)';

% --- Compute Initial Volumes and Detection Thresholds ---

[T0, T0_thT2] = Tumor_volume_th(C, L, PTV, U_0, th_T2);
[Q0, Q0_thT2] = Tumor_volume_th(C, L, PTV, Q_evo, th_T2);