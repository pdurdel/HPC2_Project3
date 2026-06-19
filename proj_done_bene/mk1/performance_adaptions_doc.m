%{

making the code mor effiecient:

(1) 'artificial_point_grad.m':
- instead of calculating every point on its own in a for loop, every valid
  points gets now calculated simultaneously

(2)  'artificial_point_grad.m':
- for every point, a small lgs gets solved, but the gradient is in a linear
  FEM-triangle constant
- now: the gradient is computed only once per trangle and then assigned to
  all evaluation points in that trangle

(3) 'artificial_points_rand.m'
- originally, each random point required a separate search for the 
  corresponding triangle
- now: all triangle assignments are computed simultaneously using cumulative 
  triangle areas and discretize
%}