function [S,Rho,M]=SepBound7StateETFSeeSaw_Random_start(NStarts)
% This function optimeze the self testing witness for 7 states ETF in dim
% 4 and the optimization is performed on separable preparation.
% Thus giving the separable preparation bound for the 7 ETF self testing
% witness.

% INPUT:
%   NStarts : number of random initializations.
%
% OUTPUT:
%   S         : largest witness value found
%   Rho       : optimal 7 product density matrices
%   M         : corresponding 21 optimal binary projective measurements

if nargin < 1
    NStarts = 10;
end

format long;
S=0;

fprintf('--------------------------------------------------\n');
    fprintf('7-state ETF separable optimization using seesaw with random input\n');
    fprintf('Number of random starts = %d\n',NStarts);
    fprintf('--------------------------------------------------\n');

for Count=1:NStarts
    % preparing 7 density matrices, as per "Appendix C 2011 New J. Phys. 13 053047"
        Rho = cell(7,1); % for storing resultant states
        Rho_sdp = cell(7,1); % for storing sdp variables for states
        Fs = [];
        for a0 = 1:7
            Rho_sdp{a0} = sdpvar(4,4,'hermitian','complex'); 
            Fs =  [Fs; Rho_sdp{a0} >= 0; trace(Rho_sdp{a0}) == 1]; % positivity and unity trace
            Fs=[Fs;PartialTranspose(Rho_sdp{a0})>=0];
            Rho{a0} = RandomDensityMatrix(4);
        end
    
    
       
    % preparing 21 binary outcome measurements.
    % specified by 2 index (i,j) with i>j and i,j={1,2,...,7}.
        M = cell(21,2); % for storing resultant measurements
        M_sdp = cell(21,2); % for storing sdp variables for measurements
        Fm = [];
        for i = 1:21
            sum = 0;
                R = RandomPOVM(4,2); % random povm initialization for the first iteration
                for k = 1:2
                    M_sdp{i,k} = sdpvar(4,4,'hermitian','complex');
                    Fm = [Fm; M_sdp{i,k} >= 0;]; % positivity
                    sum = sum + M_sdp{i,k};
                    M{i,k} = R{k}; 
                end
                Fm = [Fm; sum == eye(4);]; % completeness
        end
        vstepM = 100; vstepS = 0;Iteration=0; % declaring loop parameters 
        while (abs(vstepM-vstepS) >= 0.000000001 && Iteration<40000)
            Iteration=Iteration+1;
            
            % sdp for measurements.
            stepM = EntanglementSelfTestValue(Rho,M_sdp); % the success metric for 4 AES
           
    
            % the measurement optimization step
            diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'mosek','verbose',0));
            % preparing states for the state sdp.
            for i = 1:21
                for k = 1:2
                    M{i,k} = value(M_sdp{i,k});
                end  
            end
             % extracting the value of the objective function for the measurement step
            vstepM = value(stepM);
    
            % sdp for states
            stepS = EntanglementSelfTestValue(Rho_sdp,M); % the success metric
            % the state optimization step
            diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'mosek','verbose',0));
            % preparing states for the measurement sdp.
            for a0 = 1:7
                Rho{a0} = value(Rho_sdp{a0});
            end
            % extracting the value of the objective function for the state step
            vstepS = value(stepS);
            
        end
    S=max(S,vstepS);
    fprintf('Count:%d and Obj Value %d and the overall optimum value is %d\n ',Count,vstepS,S)
end

 fprintf('\n===============================================\n');
 fprintf('Best objective value = %.15f\n',S);
 fprintf('===============================================\n');

end

function S = EntanglementSelfTestValue(Rho,M)
% for calculating the success metric given states and measurements
    S=0;
    for i = 1:7
        for j = 1:i-1
            x=(i-1)*(i-2)/2+j;
            
            S = S + (1/(7*sqrt(14)))*(trace(Rho{i}*M{x,2})-trace(Rho{j}*M{x,2}));
        end
    end
   S=real(S); 
end

