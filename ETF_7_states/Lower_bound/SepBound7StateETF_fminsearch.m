function [bestVal, Rho, bestParam] = SepBound7StateETF_fminsearch(nStarts)
% SepBound7StateETF_fminsearch
%
% Maximizes the separable value of the 7-state ETF witness using fminsearch.
%
% The preparations are restricted, without loss of optimality, to
% pure product two-qubit states |psi_i> = |a_i> \otimes |b_i>.
% The binary measurements are optimized analytically. Therefore fminsearch 
% only optimizes the preparation states.


% INPUT:
%   nStarts : number of random initializations.
%
% OUTPUT:
%   bestVal   : largest witness value found
%   Rho       : optimal 7 product density matrices
%   M         : corresponding 21 optimal binary projective measurements
%               M{x,1}, M{x,2}
%   bestParam : optimal fminsearch parameter vector

    if nargin < 1
        nStarts = 50;
    end

    format long;

    % -------------------------------------------------------------
    % Because a common local unitary does not change the objective,
    % fix preparation 1 to |00>.
    %
    % Remaining 6 preparations: each has four parameters (thetaA, phiA,
    % thetaB, phiB). Total number of parameters = 6*4 = 24.
    % -------------------------------------------------------------

    nPar = 24;

    bestVal   = -Inf;
    bestParam = [];
    opts = optimset('Display', 'off', 'MaxIter', 200000, 'MaxFunEvals',...
        1000000, 'TolX', 1e-12, 'TolFun', 1e-14);

    fprintf('--------------------------------------------------\n');
    fprintf('7-state ETF separable optimization using fminsearch\n');
    fprintf('Number of random starts = %d\n',nStarts);
    fprintf('--------------------------------------------------\n');

    for count = 1:nStarts
        % Random initialization
        p0 = randomInitialParameters();

        % fminsearch performs minimization, therefore minimize -W
        [pOpt,negVal] = fminsearch(@objective,p0,opts);
        value = -negVal;
        if value > bestVal
            bestVal   = value;
            bestParam = pOpt;

            fprintf('Start %4d : NEW BEST = %.12f\n', ...
                    count,bestVal);
        else
            fprintf('Start %4d : %.12f   best = %.12f\n', ...
                    count,value,bestVal);
        end
    end

    % -------------------------------------------------------------
    % Construct the optimal density matrices
    % -------------------------------------------------------------

    Psi = parametersToStates(bestParam);

    Rho = cell(7,1);

    for i = 1:7
        Rho{i} = Psi{i}*Psi{i}';
    end


    

    fprintf('\n===============================================\n');
    fprintf('Best objective value = %.15f\n',bestVal);
    fprintf('===============================================\n');


    % =============================================================
    % Nested objective
    % =============================================================
    function f = objective(p)

        psi = parametersToStates(p);
        W = 0;
        coeff = 1/(7*sqrt(14));

        for ii = 1:7
            for jj = 1:ii-1

                overlap = psi{ii}'*psi{jj};
                overlapSq = abs(overlap)^2;
                traceDistance = sqrt(max(0,1-overlapSq));

                W = W + coeff*traceDistance;
            end
        end

        % fminsearch minimizes
        f = -W;
    end


    % =============================================================
    % Convert the 24 parameters into seven two-qubit product states
    % =============================================================
    function psi = parametersToStates(p)

        psi = cell(7,1);

        % Gauge fixing
        ket0 = [1;0];
        psi{1} = kron(ket0,ket0);

        ind = 1;
        for s = 2:7

            thetaA = p(ind);
            phiA   = p(ind+1);
            thetaB = p(ind+2);
            phiB   = p(ind+3);
            ind = ind + 4;

            ketA = qubitState(thetaA,phiA);
            ketB = qubitState(thetaB,phiB);

            psi{s} = kron(ketA,ketB);
        end
    end


    % =============================================================
    % Pure qubit state
    % =============================================================
    function q = qubitState(theta,phi)
        q = [cos(theta/2);
            exp(1i*phi)*sin(theta/2)];
        q = q/norm(q);
    end

    % =============================================================
    % Random initialization
    % =============================================================
    function p0 = randomInitialParameters()
        p0 = zeros(nPar,1);
        ind = 1;
        for s = 2:7
            % z = 2*rand - 1;
            % thetaA = acos(z);
            thetaA = pi*rand;
            phiA   = 2*pi*rand;

            % z = 2*rand - 1;
            % thetaB = acos(z);
            thetaB = pi*rand;
            phiB   = 2*pi*rand;

            p0(ind:ind+3) = [thetaA;phiA;thetaB;phiB];
            ind = ind + 4;
        end
    end

end
