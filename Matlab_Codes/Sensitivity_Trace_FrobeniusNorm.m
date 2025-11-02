% ========================================================================
%  Title: Cracking Composite Design: Trace-Driven Scaling for 
%         Cluster-Based Universal Design
%  Author: Francesco Danzi
%  Reference:
%     F. Danzi, "Cracking Composite Design: Trace-Driven Scaling for 
%     Cluster-Based Universal Design," Composite Structures, 2025.
%     DOI: [insert DOI once published]
%
% ------------------------------------------------------------------------
%  Description:
%     This script computes the sensitivity as detailed in Section 3.2 
%     of the manuscript. In particular, it evaluates the sensitivity 
%     with respect to the Frobenius Norm.
%
%  License:
%     MIT License (recommended for research sharing)
%
%  GitHub Repository:
%     https://github.com/FDanzi/Clustering-Composite-Materials
% ========================================================================


clc; clear; close all;

%% Define Fiber Properties [E1 (GPa), E2 (GPa), G12 (GPa), nu12, Density]
fibers = {
    'High-Modulus Carbon', [540, 15, 23, 0.2, 1.8];
    'Intermediate-Modulus Carbon', [294, 10, 14, 0.3, 1.8];
    'Standard-Modulus Carbon', [230, 15, 25, 0.2, 1.76];
    'E-Glass', [72, 72, 30, 0.22, 2.54];
    'S-Glass', [86, 86, 35, 0.23, 2.49];
    'Boron', [400, 40, 18, 0.2, 2.6];
    'Kevlar 49', [130, 5, 2.5, 0.34, 1.44];
    'Kevlar 29', [80, 4, 2.0, 0.36, 1.44];
     % Natural Fibers
    'Flax',[60,   12,   4,   0.3,  1.4];   % Flax (using mid-range values)
    'Hemp',[45,   10,   3.25, 0.32, 1.45]; % Hemp
    'Jute',[35,    7.5, 2.5, 0.35, 1.3];   % Jute
    'Kenaf',[42,    8.5, 3.1, 0.34, 1.4];   % Kenaf
    'Ramie',[55,   12.5, 4,   0.3,  1.5];   % Ramie
    'Coir (Coconut)',[5,     1.5, 0.75, 0.4,  1.2];  % Coir (Coconut)
};

%% Define Matrix Properties [E (GPa), G (GPa), nu, Density]
matrices = {
    'Epoxy', [3.5, 1.3, 0.35, 1.2];
    'Bismaleimide (BMI)', [3.8, 1.4, 0.35, 1.25];
    'Polyimide', [2.8, 1.0, 0.34, 1.4];
    'PEEK', [3.7, 1.3, 0.38, 1.3];
    'PEI', [3.2, 1.2, 0.38, 1.27];
    'PPS', [3.5, 1.3, 0.37, 1.35];
};

%% Fiber Volume Fractions
Vf_range = 0.4:0.0375:0.7; % 40% to 70%

%% Small perturbation for finite difference
dPi = 1e-4;

%% Function to compute Q* matrix
compute_Qstar = @(p1, p2, p3) (1 / (1 + p1 + 2*p2*(1 - p3^2*p1))) * ...
    [1, p2, p1;
     p2, p1, 0;
     0, 0, 2*(1 - p3^2*p1)*p2];

%% Initialize Results Storage
results = [];

%% Loop Over All Fiber-Matrix Combinations
for f = 1:size(fibers,1)
    fiber_name = fibers{f,1};
    Ef1 = fibers{f,2}(1);
    Ef2 = fibers{f,2}(2);
    Gf12 = fibers{f,2}(3);
    nuf12 = fibers{f,2}(4);

    for m = 1:size(matrices,1)
        matrix_name = matrices{m,1};
        Em = matrices{m,2}(1);
        Gm = matrices{m,2}(2);
        num = matrices{m,2}(3);

        for Vf = Vf_range
            Vm = 1 - Vf;

            %% Compute Homogenized Composite Properties
            E1 = Vf * Ef1 + Vm * Em; % Rule of Mixtures
            etaE=(Ef2/Em-1)/(Ef2/Em+2);
            E2 = Em * (1+2*etaE*Vf)/(1-etaE*Vf); % Halpin-Tsai
            etaG=(Gf12/Gm-1)/(Gf12/Gm+2);
            G12 = Gm * (1 + 1*etaG*Vf) / (1 - etaG*Vf); % Halpin-Tsai
            nu12 = Vf * nuf12 + Vm * num; % Rule of Mixtures

            %% Compute Non-Dimensional Parameters
            pi1 = E2 / E1;
            pi2 = G12 / E1;
            pi3 = nu12;

            %% Compute Baseline Q* Matrix
            Q_base = compute_Qstar(pi1, pi2, pi3);

            %% Compute Sensitivities
            dF_pi1 = norm(compute_Qstar(pi1 + dPi, pi2, pi3) - Q_base, 'fro') / dPi;
            dF_pi2 = norm(compute_Qstar(pi1, pi2 + dPi, pi3) - Q_base, 'fro') / dPi;
            dF_pi3 = norm(compute_Qstar(pi1, pi2, pi3 + dPi) - Q_base, 'fro') / dPi;

            %% Identify Most Influential π-Group
            [~, maxIndex] = max([dF_pi1, dF_pi2, dF_pi3]);
            pi_labels = ["E2/E1", "G12/E1", "ν12"];
            dominant_pi = pi_labels(maxIndex);

            %% Store Results
            results = [results; {fiber_name, matrix_name, Vf, pi1, pi2, pi3, dF_pi1, dF_pi2, dF_pi3, dominant_pi}];
        end
    end
end

%% Convert Results to Table
results_table = cell2table(results, ...
    'VariableNames', {'Fiber', 'Matrix', 'Vf', 'Pi1 (E2/E1)', 'Pi2 (G12/E1)', 'Pi3 (ν12)', ...
                      'Sensitivity Pi1', 'Sensitivity Pi2', 'Sensitivity Pi3', 'Dominant Pi'});

%% Display Results Table
disp(results_table);

%% Save Table to CSV
writetable(results_table, 'sensitivity_analysis_results.csv');
