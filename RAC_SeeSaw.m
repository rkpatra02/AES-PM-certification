function [Initial_Value,Rho,M]=RAC_SeeSaw()
% This function calculate the RAC success probability with states from
% C2*C2. 
% The input on Alice side are two dits x=x_0x_1={0,1,2,3}^{*2}
% Bob has two input y={0,1} and produce an outcome dit b={0,1,2,3}
% Its perform huristic see saw search over separable encoding.

format long;

%--------------Initial parameter------------------------------

Nd=4; % Dimension of the input string
ND=4; % Dimension of Alice to Bob message 

%-------------Number of initial randomization-----------
N_rand=10;
Initial_Value=0;
for m=1:N_rand
    % ----------preparing Nx no density matrices, as preparation---------------
        Nx=Nd^2; % Number of preparation
        Rho = cell(Nx,1); % For storing resultant states during Seesaw
        RhoFinal = cell(Nx,1); % For storing the final optimal states after all random selection
        Rho_sdp = cell(Nx,1); % For storing sdp variables for states
        Fs = [];
        for a0 = 1:Nx
            Rho_sdp{a0} = sdpvar(ND,ND,'hermitian','complex'); 
            Fs =  [Fs; Rho_sdp{a0} >= 0; trace(Rho_sdp{a0}) == 1]; % Positivity and unity trace
            Fs=[Fs;PartialTranspose(Rho_sdp{a0})>=0]; % Separability condition imposed by PPT
            Rho{a0} = RandomDensityMatrix(ND,1); % Storing random density matrix for See-saw iteration
        end
    
       
    % -------------preparing 'Nm' measurement with 'No' outcomes.-------------
    % 
        Nm=2; % Number of measurement
        No=Nd; % Number of outcome of each measurement
        M = cell(Nm,No); % For storing resultant measurements during seesaw
        MFinal = cell(Nm,No); % For storing final resultant measurements over all initial random selection
        M_sdp = cell(Nm,No); % For storing sdp variables for measurements
        Fm = [];
        for i = 1:Nm
            sum = 0;
                R = RandomPOVM(ND,No); % Random povm initialization for the first iteration
                for k = 1:No
                    M_sdp{i,k} = sdpvar(ND,ND,'hermitian','complex');
                    Fm = [Fm; M_sdp{i,k} >= 0;]; % Positivity
                    sum = sum + M_sdp{i,k};
                    M{i,k} = R{k}; 
                end
                Fm = [Fm; sum == eye(ND);]; % Completeness of measurement
        end
    
        %---------costruction of Seesaw loop-------------------------------
    
        
    
        vstepM = 100; vstepS = 0;Iteration=0; % Declaring loop parameters 
        
        while (abs(vstepM-vstepS) >= 0.000000001 && Iteration<1000)
            Iteration=Iteration+1;
            
            % SDP optimization of measurements.
            stepM = RacSuccess(Rho,M_sdp); % the success metric for 4 AES
           
            diagnostics = optimize([Fm;], -stepM, sdpsettings('solver', 'mosek','verbose',0));
            % Preparing Measurements for the state sdp.
            for i = 1:Nm
                for k = 1:No
                    M{i,k} = value(M_sdp{i,k});
                end  
            end
             % Extracting the value of the objective function for the measurement step
            vstepM = value(stepM);
    
            % SDP optimization of states
            stepS = RacSuccess(Rho_sdp,M); % the success metric
            
            diagnostics = optimize([Fs;], -stepS, sdpsettings('solver', 'mosek','verbose',0));
            % Preparing states for the measurement sdp.
            for a0 = 1:Nx
                Rho{a0} = value(Rho_sdp{a0});
            end
            % Extracting the value of the objective function for the state step
            vstepS = value(stepS);
        end

        if vstepS>Initial_Value
            Initial_Value=vstepS;
            RhoFinal=Rho;
            MFinal=M;
        end
        fprintf('Optimal value for randomization %d is %d \n',m,vstepS)
end
fprintf('Final optimal value with %d randomization is %d\n',N_rand,Initial_Value)
end


function S = RacSuccess(Rho,M)
% for calculating the success metric given states and measurements
    
    S=0;Nx=4;Ny=2;Nb=4;
    for x1 = 1:Nx
        for x2 = 1:Nx
            % size(Rho{Nx*(x1-1)+x2})
            % size(M{1,x1})
            S=S+trace(Rho{Nx*(x1-1)+x2}*M{1,x1});
            S=S+trace(Rho{Nx*(x1-1)+x2}*M{2,x2});
        end
    end
   S=S/(Nx^2*Ny); 
end
