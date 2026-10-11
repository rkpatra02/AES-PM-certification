function [S,Rho,M]=SepBound16StateSICSeeSaw_Warmstart(Nstarts)

% SepBound16StateSICSeeSaw_Warmstart
%
% This function computes a lower bound on the separable value of the
% self-testing witness associated with 16 SIC states in dimension 4.
% The optimization is performed over separable two-qubit preparations
% and binary measurements using an alternating see-saw procedure.
%
% The preparation step is formulated as an SDP using the PPT criterion,
% which is necessary and sufficient for separability in the 2 x 2 case.
%
% The see-saw is warm-started from a set of near-optimal separable
% preparations stored in the file OptSIC_Rho.mat.
%
% -------------------------------------------------------------------------
% INPUT:
%   Nstarts  - Number of independent see-saw runs.
%
% OUTPUT:
%   S        - Largest objective value obtained over all runs.
%   Rho      - Preparation states from the final run.
%   M        - Measurement operators from the final run.
% -------------------------------------------------------------------------

format long;
S=0;
% Load the warm-start preparation states stored in OptSIC_Rho.mat.
data = load('OptSIC_Rho.mat');
OptSIC_Rho = data.OptSIC_Rho;
if nargin<1
    Nstarts=2;
end

fprintf('--------------------------------------------------\n');
fprintf('16-state SIC witness: separable optimization using warm-started see-saw\n');
fprintf('Number of see-saw runs = %d\n',Nstarts);
fprintf('--------------------------------------------------\n');



for Count=1:Nstarts
        
    % ---------------------------------------------------------------------
    % Define the 16 separable preparation states.
    %
    % Rho stores the current numerical preparation states.
    % Rho_sdp stores the corresponding SDP variables.
    % In the two-qubit case, positivity, unit trace, and positivity under
    % partial transpose characterize the separable-state set exactly.
    % ---------------------------------------------------------------------
    
        Rho = cell(16,1); % For storing resultant states
        Rho_sdp = cell(16,1); % For storing sdp variables for states
        Fs = [];
        for a0 = 1:16
            Rho_sdp{a0} = sdpvar(4,4,'hermitian','complex'); 
            % Positivity and unit-trace constraints.
            Fs =  [Fs; Rho_sdp{a0} >= 0; trace(Rho_sdp{a0}) == 1]; 
            % PPT constraint; exact for separability in dimension 2 x 2.
            Fs=[Fs;PartialTranspose(Rho_sdp{a0})>=0];
            % Warm-start from the stored near-optimal preparation states.
            Rho{a0} = OptSIC_Rho(:,:,a0);
        end
    
    
       
    % ---------------------------------------------------------------------
    % Define the 120 binary measurements.
    %
    % Each measurement is associated with an unordered pair (i,j),
    % with i > j and i,j in {1,...,16}. Since C(16,2) = 120, there are
    % 120 measurement settings, each with two outcomes.
    % ---------------------------------------------------------------------
        M = cell(120,2); % For storing resultant measurements
        M_sdp = cell(120,2); % For storing sdp variables of measurements
        Fm = [];
        for i = 1:120
            sum = 0;
                R = RandomPOVM(4,2); % Random POVM initialization for the first step.
                for k = 1:2
                    M_sdp{i,k} = sdpvar(4,4,'hermitian','complex');
                    Fm = [Fm; M_sdp{i,k} >= 0;]; % % Positivity of each POVM element.
                    sum = sum + M_sdp{i,k};
                    M{i,k} = R{k}; 
                end
                Fm = [Fm; sum == eye(4);]; % % Completeness of the POVM.
        end

    % ---------------------------------------------------------------------
    % Alternate between measurement optimization and preparation
    % optimization until convergence or until the maximum number of
    % iterations is reached.
    % ---------------------------------------------------------------------

        % Declaring loop parameters 
        vstepM = 100; 
        vstepS = 0;
        Iteration=0; 

        while (abs(vstepM-vstepS) >= 0.000000001 && Iteration<2000)
            Iteration=Iteration+1;
            
        % -----------------------------------------------------------------
        % Measurement optimization step:
        % optimize all binary POVMs while keeping the preparations fixed.
        % -----------------------------------------------------------------
            
            stepM = EntanglementSelfTestValue(Rho,M_sdp); % the success metric for 4 AES
           
            diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'mosek','verbose',0));

            % Extract the optimized measurement operators.
            for i = 1:120
                for k = 1:2
                    M{i,k} = value(M_sdp{i,k});
                end  
            end

            % Objective value after the measurement step.
            vstepM = value(stepM);
    
        % -----------------------------------------------------------------
        % Preparation optimization step:
        % optimize all separable preparation states while keeping the
        % measurements fixed.
        % -----------------------------------------------------------------

            stepS = EntanglementSelfTestValue(Rho_sdp,M); % the success metric
            diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'mosek','verbose',0));

            % Extract the optimized preparation states
            for a0 = 1:16
                Rho{a0} = value(Rho_sdp{a0});
            end
            % Objective value after the preparation step.
            vstepS = value(stepS);
            
        end
    S=max(S,vstepS);
    fprintf('Run %d: objective value = %.12f, best value = %.12f\n', ...
            Count,vstepS,S);
end

fprintf('\n===============================================\n');
 fprintf('Best objective value = %.15f\n',S);
 fprintf('===============================================\n');
end

function S = EntanglementSelfTestValue(Rho,M)
% -------------------------------------------------------------------------
% Compute the 16-state SIC self-testing witness value for a fixed set of
% preparation states and binary measurements.
%
% The coefficient is
%
%   alpha_i alpha_j ||P_i-P_j||_1
%   = (1/16)^2 * 4/sqrt(5),
% -------------------------------------------------------------------------
    S=0;
    for i = 1:16
        for j = 1:i-1
            x=(i-1)*(i-2)/2+j;
            S = S + (1/16)^2*(4/sqrt(5))*(trace(Rho{i}*M{x,2})-trace(Rho{j}*M{x,2}));
        end
    end
  end

