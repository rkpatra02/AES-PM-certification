% This script calculates an upper bound on the separable value of the RAC
% witness using an SDP relaxation with projective measurements.
% The script uses the problem definition contained in RAC44sepUB_Proj.m.

settings = NVSettings;
problem = RAC44sepUB_Proj;

fprintf('---------------------------------------------------------------------------\n')
fprintf('Upper bound on the separable value of the RAC witness using an SDP relaxation with projective measurement\n')
fprintf('---------------------------------------------------------------------------\n')
% ------------------------------------------------------------------------------
% Different monomial sets correspond to different relaxation levels.
% Choose the desired level by selecting the appropriate monomial set below.
% ------------------------------------------------------------------------------
monomials = {'families' [] [1] [2] [1 2] [2 1] [2 2]};
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2]};
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2] [2 2 2]};

format long

%nvOptimize(problem, monomials, 'none', settings)
%nvOptimize(problem, monomials, 'reynolds', settings)
%nvOptimize(problem, monomials, 'isotypic', settings)
%nvOptimize(problem, monomials, 'irreps', settings)
nvOptimize(problem, monomials, 'blocks', settings)
