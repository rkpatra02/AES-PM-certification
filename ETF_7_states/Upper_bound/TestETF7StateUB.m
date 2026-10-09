

settings = NVSettings;
problem = ETF7StateUB;
%-------------------------------------------------------------------------------
%    Different monomials corresponds to different levels choose accordingly
% ------------------------------------------------------------------------------
 monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2]};
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2] [1 1 1]};
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2] [1 1 1] [1 1 2] };
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2] [1 1 1] [1 1 2] [2 1 1]};
%monomials = {'families' [] [1] [2] [1 1] [1 2] [2 1] [2 2] [1 1 1] [1 1 2] [1 2 2]};

format long
%nvOptimize(problem, monomials, 'none', settings)
%nvOptimize(problem, monomials, 'reynolds', settings)
%nvOptimize(problem, monomials, 'isotypic', settings)
%nvOptimize(problem, monomials, 'irreps', settings)
nvOptimize(problem, monomials, 'blocks', settings)


% SEP THRESHOLD: Level 2 -> 0.75


