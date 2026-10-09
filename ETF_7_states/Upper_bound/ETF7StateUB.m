classdef ETF7StateUB < NVProblem
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self)
        % The first 7 operators are the states, the next 7 operators are 
        % the partial A states, the next 7 operators are the partial B states, 
        % the next 7*6/2 = 21 are the traceless observables (two-outcome measurements)
            dim = 4;
            % Total operators -> 21 + 21 = 42
            X = cell(1, 42);
     	    XA = cell(1, 7);
	        XB = cell(1, 7);
                for i = 1:7
                    XA{i} = qdimsum.Random.pureNormalizedDensityMatrix(2);
                    XB{i} = qdimsum.Random.pureNormalizedDensityMatrix(2);
                    X{i} = kron(XA{i},XB{i});
                end
	        for i = 1:7
	        X{i+7} = kron(XA{i},eye(2));
	        X{i+14} = kron(eye(2),XB{i}); 
            end
            y = 21;
            for x1 = 1:7
                for x2 = x1+1:7
                U = qdimsum.Random.unitary(4);
                y = y + 1;
                X{y} = U*[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]*U';
                end;
            end
        end
        function C = operatorSDPConstraints(self, X)
            id = eye(4);
            C = cell(1, 42);
            for i = 1:21
            C{2*i-1} = id + X{21+i};   
            C{2*i  } = id - X{21+i};   
            end
        end
        function K = sampleStateKraus(self)
            K = eye(4); % dimension is 4
        end
        function obj = computeObjective(self, X, K)
            obj = 0;
            y = 22;
            for x1 = 1:7
                for x2 = x1+1:7
                        rho = X{x1};
                        M = X{y};
                        obj = obj + real(trace(M * rho));
                        rho = X{x2};
                        obj = obj - real(trace(M * rho));
                        y = y + 1;
                end
            end
            obj = obj/2;
            obj = (1/(7*sqrt(14)))*obj;
        end
        function generators = symmetryGroupGenerators(self)
        % Minimal generators (length 42) for V7sepcorr that preserve the objective.
        
        % p1 : transposition 1 <-> 2 (sign is on #22)
p1 = [ ...
  2,  1,  3,  4,  5,  6,  7, ...
  9,  8, 10, 11, 12, 13, 14, ...
 16, 15, 17, 18, 19, 20, 21, ...
 -22, 28, 29, 30, 31, 32, 23, ...
 24, 25, 26, 27, 33, 34, 35, 36, ...
 37, 38, 39, 40, 41, 42 ];

% pcycle : 7-cycle 1->2->3->4->5->6->7->1 (signs are on #22, #23, #24, #25, #26 and #27)
pcycle = [ ...
  2,  3,  4,  5,  6,  7,  1, ...
  9, 10, 11, 12, 13, 14,  8, ...
 16, 17, 18, 19, 20, 21, 15, ...
 28, 29, 30, 31, 32, -22, 33, ...
 34, 35, 36, -23, 37, 38, 39, ...
 -24, 40, 41, -25, 42, -26, -27 ];
        
        generators = [p1; pcycle];

    end; 
        
        function types = operatorTypes(self)
            types = {1:21 22:42};
        end
    end
        end