function []=RAC_Parametric_plot()
% -------------------------------------------------------------------------
% Construct the 16 entangled preparation states.
% The states are generated from a fiducial state using the generalized
% four-dimensional shift and phase operators PX and PZ.
% -------------------------------------------------------------------------
psi=[sqrt(3)/2,1/(2*sqrt(3)),1/(2*sqrt(3)),1/(2*sqrt(3))]';
PX=[0,0,0,1;1,0,0,0;0,1,0,0;0,0,1,0];     % Generalized shift operator in dimension 4
PZ=[1,0,0,0;0,1j,0,0;0,0,-1,0;0,0,0,-1j]; % Generalized phase operator in dimension 4
for x0=0:3
    for x1=0:3
        psi_state_ENT{4*x0+x1+1}=PX^(x0)*PZ^(x1)*psi;
    end
end

% -------------------------------------------------------------------------
% Construct the corresponding 16 separable preparation states.
% -------------------------------------------------------------------------
a=cos(pi/8);b=sin(pi/8);
X = [0 1; 1 0];
Z = [1 0; 0 -1];
P = [1 0; 0 1i];

chi=kron([a;b],[a;b]);
F=[1 1 1 1; 1 1 1i 1; 1i -1i 1 -1; 1i -1 1i 1]; % Phase factors used in the separable-state construction.
M=[0 0;0 1;1 0;1 1]; % Binary representation of x0 used to generate the local X operations.

% -------------------------------------------------------------------------
% Define the range of the interpolation parameter c.
% The parametric preparation states are
%
% |psi_x(c)> = c |psi_x^ENT>
%            + sqrt(1-c^2) |psi_x^SEP>,
%
% followed by normalization.
% -------------------------------------------------------------------------

for x0=0:3
    for x1=0:3
        psi_state_sep{4*x0+x1+1}=F(x0+1,x1+1)*(kron(Z,P)^x1)*(kron(X^(M(x0+1,1)),X^(M(x0+1,2))))*chi;
    end
end

    N=200; % Number of intervals
    x1=0;
    x2=1;
    Yval=zeros(N+1,1);
    Xval=x1:(x2-x1)/N:x2;

% -------------------------------------------------------------------------
% For each value of c, construct the 16 parametric preparation states and
% optimize Bob's measurements to obtain the corresponding RAC witness value.
% -------------------------------------------------------------------------
    for i=1:N+1
        c=Xval(i);
         for j=1:16
            psi_state_final{j}=c*psi_state_ENT{j}+sqrt(1-c^2)*psi_state_sep{j};
            % Normalize the interpolated state.
            psi_state_final{j}=psi_state_final{j}/norm(psi_state_final{j});
            % Convert the state vector to a density operator.
            Psi_density_final{j}=psi_state_final{j}*psi_state_final{j}';
        end
        % Optimize the RAC witness violation for the current preparation ensemble.
        Yval(i)=OptimalMeasurementSDP(Psi_density_final);
        
    end

% -------------------------------------------------------------------------
% Separable RAC bound.
% -------------------------------------------------------------------------
  
pSuc = (1/4)*(1+1/sqrt(2))^2;
m = pSuc * ones(size(Xval));

% Convert arrays to column vectors for plotting.
Xval = Xval(:);
Yval = Yval(:);
m    = m(:);



% -------------------------------------------------------------------------
% Plot the RAC witness value as a function of c and indicate the region in
% which the separable bound is violated.
% -------------------------------------------------------------------------
figure;

plot(Xval, Yval, 'LineWidth', 2);
hold on;

plot(Xval, m, '--', 'LineWidth', 2);

% Identify the parameter region where W^RAC > Q^SEP.
idx = Yval > m;

xpatch = [Xval(idx); flipud(Xval(idx))];
ypatch = [Yval(idx); flipud(m(idx))];

fill(xpatch, ypatch, 'green', ...
    'FaceAlpha', 0.3, ...
    'EdgeColor', 'none');

% labels (only these fonts increased)
xlabel('c', ...
       'Interpreter','latex', ...
       'FontSize',18);

ylabel('$W^{\mathrm{RAC}}$', ...
       'Interpreter','latex', ...
       'FontSize',18);

% legend (font increased only here)
legend({'$Q$', '$Q^{\mathrm{SEP}}$', '$Q > Q^{\mathrm{SEP}}$'}, ...
       'Interpreter','latex', ...
       'FontSize',16, ...
       'Location','best');

grid on;
hold off;

end



function M=OptimalMeasurementSDP(Rho)
% -------------------------------------------------------------------------
% Optimize Bob's two four-outcome POVMs for a fixed set of preparation
% states by solving a semidefinite program.
% This function return the optimal RAC witness violation for the fixed
% preparations.
% -------------------------------------------------------------------------
    M_sdp = cell(2,4); % SDP variables for the measurement operators
    Fm = [];
    for i = 1:2
        sum = 0;
           
            for k = 1:4
                M_sdp{i,k} = sdpvar(4,4,'hermitian','complex');
                % Positivity of each POVM element.
                Fm = [Fm; M_sdp{i,k} >= 0;]; 
                sum = sum + M_sdp{i,k};
                 
            end
            % Completeness of each POVM.
            Fm = [Fm; sum == eye(4);]; 
    end
    stepM = RacSuccess(Rho,M_sdp);
    diagonistic=optimize([Fm;], -stepM, sdpsettings('solver', 'mosek','verbose',0));
    M=value(stepM);

end
function S = RacSuccess(Rho,M)
% -------------------------------------------------------------------------
% Compute the average success probability of the 2-to-1 random access code
% with alphabet size 4.
% -------------------------------------------------------------------------
    
    S=0;Nx=4;Ny=2;Nb=4;
    for x1 = 1:Nx
        for x2 = 1:Nx
            S=S+trace(Rho{Nx*(x1-1)+x2}*M{1,x1});
            S=S+trace(Rho{Nx*(x1-1)+x2}*M{2,x2});
        end
    end
   S=S/(Nx^2*Ny); 
   S=real(S);
end
