function [S,RhoOpt,MOpt]=SepBound16StateSICSeeSaw_Randomstart(Nstarts)
% This function optimeze the self testing witness for 16 SIC states in dim
% 4 and the optimization is performed on separable preparation.
% Thus giving the separable preparation bound for the 16 SIC self testing
% witness.
% This randomized the intial input Nstarts times and return the overall
% optimal value.
%---------------------------------------
% Input:
% Nstarts = Number of random intializatiom
%
% Output:
% S = Optimal objective value
% RhoOpt = Optimal separable preparations
% MOpt = Optimal measurements
%-----------------------------

format long;
S=0;
if nargin<1
    Nstarts=10;
end

fprintf('--------------------------------------------------\n');
fprintf('16-state SIC witness: separable optimization using randomized input see-saw\n');
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
        Rho = cell(16,1); % for storing resultant states
        Rho_sdp = cell(16,1); % for storing sdp variables for states
        Fs = [];
        for a0 = 1:16
            Rho_sdp{a0} = sdpvar(4,4,'hermitian','complex'); 
            % Positivity and unity trace
            Fs =  [Fs; Rho_sdp{a0} >= 0; trace(Rho_sdp{a0}) == 1]; 
            % PPT constraint; exact for separability in dimension 2 x 2.
            Fs=[Fs;PartialTranspose(Rho_sdp{a0})>=0];
            % Random preparation states intialization
            Rho{a0} = RandomDensityMatrix(4,1);
        end
    
    
       
    % ---------------------------------------------------------------------
    % Define the 120 binary measurements.
    %
    % Each measurement is associated with an unordered pair (i,j),
    % with i > j and i,j in {1,...,16}. Since C(16,2) = 120, there are
    % 120 measurement settings, each with two outcomes.
    % ---------------------------------------------------------------------
        M = cell(120,2); % For storing resultant measurements
        M_sdp = cell(120,2); % For storing sdp variables for measurements
        Fm = [];
        for i = 1:120
            sum = 0;
                R = RandomPOVM(4,2); % Random POVM initialization for the first step.
                for k = 1:2
                    M_sdp{i,k} = sdpvar(4,4,'hermitian','complex');
                    Fm = [Fm; M_sdp{i,k} >= 0;]; % Positivity of each POVM element.
                    sum = sum + M_sdp{i,k};
                    M{i,k} = R{k}; 
                end
                Fm = [Fm; sum == eye(4);]; % Completeness of the POVM.
        end

        % Declaring loop parameters 
        vstepM = 100; 
        vstepS = 0;
        Iteration=0; 

        while (abs(vstepM-vstepS) >= 0.000000001 && Iteration<1000)
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
        if S>vstepS
            S=vstepS;
            RhoOpt=Rho;
            MOpt=M;
            fprintf('Run %d: New best value = %.12f\n', ...
            Count,S);
        else
            fprintf('Run %d: objective value = %.12f, best value = %.12f\n', ...
            Count,vstepS,S);
        end

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

