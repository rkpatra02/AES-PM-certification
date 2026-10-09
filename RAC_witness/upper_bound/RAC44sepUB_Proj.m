classdef RAC44sepUB_Proj < NVProblem
    properties
       forceReal = true;
    end
    methods
        function X = sampleOperators(self)
        % The first 48 operators are associated with the 16 preparations:
        %   X{1:16}   = rho_x = rho_x^A \otimes rho_x^B,
        %   X{17:32}  = rho_x^A \otimes I,
        %   X{33:48}  = I \otimes rho_x^B.
        %
        % The next four operators, X{49:52}, correspond to the four
        % projectors of the first projective measurement, and the last
        % four operators, X{53:56}, correspond to the second projective
        % measurement.
            dim = 4;
            % Total number of operators: 48 preparation-related operators
            % plus 8 measurement operators, giving 56 operators in total.
            X = cell(1, 56);
            XA = cell(1, 16);
            XB = cell(1, 16);
            for i = 1:16
                XA{i} = qdimsum.Random.pureNormalizedDensityMatrix(2);
                XB{i} = qdimsum.Random.pureNormalizedDensityMatrix(2);
                X{i} = kron(XA{i},XB{i});
            end
	        for i = 1:16
	        X{i+16} = kron(XA{i},eye(2));
	        X{i+32} = kron(eye(2),XB{i}); 
            end   

            % First four-outcome projective measurement.
            U = qdimsum.Random.unitary(4);
            X{49} = U*[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]*U';
            X{50} = U*[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]*U';
            X{51} = U*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0]*U';
            X{52} = U*[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]*U';

            % Second four-outcome projective measurement.
            U = qdimsum.Random.unitary(4);
            X{53} = U*[1 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0]*U';
            X{54} = U*[0 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0]*U';
            X{55} = U*[0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 0]*U';
            X{56} = U*[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1]*U';
        end
        function K = sampleStateKraus(self)
            K = eye(4); % The Hilbert space dimension is 4
        end
        function obj = computeObjective(self, X, K)
            obj = 0;
            for x1 = 1:4
                for x2 = 1:4
                    for y = 1:2
                        if y == 1
                            b = x1;
                        else
                            b = x2;
                        end
                        rho = X{x1+(x2-1)*4};
                        M = X{3*16+b+(y-1)*4};
                        obj = obj + real(trace(M * rho)/32);
                    end
                end
            end
        end
        function generators = symmetryGroupGenerators(self)
        % Minimal generating set for the symmetry group acting on the
        % 56 operators.

        % Swap the preparation labels x1 <-> x2. The same permutation is
        % applied to the 16 preparation states and to the corresponding
        % local operators, while the two measurement blocks
        % X{49:52} and X{53:56} are exchanged.
    
        swapX1X2 = [ ...
                    1,  5,   9, 13,  2,  6, 10, 14, ...
                    3,  7,  11, 15,  4,  8, 12, 16, ...
                    17, 21, 25, 29, 18, 22, 26, 30, ...
                    19, 23, 27, 31, 20, 24, 28, 32, ...
                    33, 37, 41, 45, 34, 38, 42, 46, ...
                    35, 39, 43, 47, 36, 40, 44, 48, ...
                    53, 54, 55, 56, 49, 50, 51, 52];

        % Apply the cyclic permutation x1: 1 -> 2 -> 3 -> 4 -> 1 to all
        % preparation states and the corresponding local operators.
        % The first measurement block, X{49:52}, is permuted in the same
        % way, while the second measurement block, X{53:56}, is unchanged.

        px1_cycle = [ ...
                     2,  3,  4,  1,  6,  7,  8,  5, ...
                    10, 11, 12,  9, 14, 15, 16, 13, ...
                    18, 19, 20, 17, 22, 23, 24, 21, ...
                    26, 27, 28, 25, 30, 31, 32, 29, ...
                    34, 35, 36, 33, 38, 39, 40, 37, ...
                    42, 43, 44, 41, 46, 47, 48, 45, ...
                    50, 51, 52, 49, 53, 54, 55, 56]; 
  
        % Swap the values x1 = 1 and x1 = 2 for every x2. The same
        % transposition is applied to the corresponding local operators
        % and to outcomes 1 and 2 of the first measurement. Operators
        % associated with x1 = 3,4 and the second measurement remain
        % unchanged.            
                
        px1_trans = [ ...
                     2,  1,  3,  4,  6,  5,  7,  8, ... 
                    10,  9, 11, 12, 14, 13, 15, 16, ...
                    18, 17, 19, 20, 22, 21, 23, 24, ...
                    26, 25, 27, 28, 30, 29, 31, 32, ...
                    34, 33, 35, 36, 38, 37, 39, 40, ...
                    42, 41, 43, 44, 46, 45, 47, 48, ...
                    50, 49, 51, 52, 53, 54, 55, 56];
  
       generators = [ swapX1X2; px1_cycle; px1_trans ];
end
        function types = operatorTypes(self)
            types = {1:48 49:56};
        end
    end
end