# Set up problems of the form:
# 
# min  c^Tx + d^Tw
# st   Ax + Bw >= b
#      0 <= w <= u
#      x in argmin_x { 0.5 x^T Qx + w^Tx : Hx >= g}.
# 
# Equivalently:
# min  c^Tx + d^Tw
# st   Ax + Bw >= b
#      0 <= w <= u
#      Qx + w - H^Ty = 0
#      0 <= y \perp Hx - g >= 0.
# 
# Construct the instances in format_3.
# 
# The matrices are stored in a sparse manner:
# [[#rows,#cols,#nonzeroes],
# [0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
# #nonzeroes in rows 1 and 2 and 3,...,
# #nonzeroes in rows 1..(#rows-1)],
# [#nonzeroes in row 1,#nonzeroes in row 2,
# #nonzeroes in row 3,...,
# #nonzeroes in row (#rows-1),#nonzeroes in last row],
# [(column of each nonzero) - 1],
# [the nonzeroes]
# ]
#
# Take A, B, H dense to simplify the calculation of
# sparsity patterns.
#
# Use the QPnonconvex/bilevel problems as an initial template.
#
# Generate everything as uniform01.
#  Include nonnegativity constraints in the inner problem.
#  Any w with Bw >= b and 0<=w<=u will then be feasible.
#
# Take u=vector of ones.
#
param n > 0;   # dimension of x and w (small)
param bdim >= 0;  # dimension of b
param gdim >= 0;  # dimension of g without the x>=0 constraints.
param rank >= 0;  # Q=LL^T and L has rank columns.

param d{j in 1..n}:=(Uniform(0,1));   #  objective function
param c{j in 1..n}:=(Uniform(0,1));   #  objective function

param b{j in 1..bdim}:=(Uniform(0,1));
param g{j in 1..gdim}:=(Uniform(0,1));

param u{j in 1..n} = 1;

param A{i in 1..bdim,j in 1..n}=(Uniform(0,1)) ; # 

param B{i in 1..bdim,j in 1..n}=(Uniform(0,1)) ; # 

param H{i in 1..gdim,j in 1..n}=(Uniform(0,1)) ; # 

param L{i in 1..n,j in 1..rank}=(Uniform(-1,1)) ; # 

param Q{i in 1..n, j in 1..n} := sum{k in 1..rank} L[i,k]*L[j,k];
#for {i in 1..n} {
#  for {j in 1..n} {
#    let Q[i,j] := sum{k in 1..rank} L[i,k]*L[j,k];
#  }
#}

#
var x{j in 1..n} >= 0;
var w{j in 1..n} >= 0, <= u[j];

var y{j in 1..gdim+n} >= 0;
var ytilde{j in 1..n};

var sigma{j in 1..n};
#

minimize bilevelobj: sum{i in 1..n} (c[i] * x[i] + d[i]*w[i]);

subject to bcon{i in 1..bdim}:
  sum{l in 1..n} A[i,l]*x[l] + sum{k in 1..n} B[i,k]*w[k]  >= b[i];

subject to gcon{i in 1..gdim}:
  sum{j in 1..n} H[i,j]*x[j] >= g[i];

subject to def_ytilde{i in 1..n}:
  ytilde[i] = y[gdim+i] + sum{j in 1..gdim} H[j,i]*y[j];

subject to KKTgrad{i in 1..n}:
  sum{j in 1..n} Q[j,i]*x[j] + w[i] - ytilde[i] = 0;

subject to lifted:
  sum{i in 1..n} sigma[i] - sum{j in 1..gdim} g[j]*y[j] <= 0;

#subject to quadcon{i in 1..n}:
#  (x[i] + ytilde[i])^2 <= 4*sigma[i] + u[i]*(ytilde[i]-x[i]);
#
# to use the quadcon, need to incorporate Q.
# Eg, if L^{-1} exists, then we have
#   w = L (L^{-1} ytilde - L^Tx)
# and we can choose
#   sigma = (L^{-1} ytilde) \bullet (L^Tx)
#
# If in addition every entry of L^{-1} is nonnegative then L^{-1}w
# is bounded above by u.
#
# Then quadcon can be expressed as
#  ((L^Tx)[i] + (L^{-1}ytilde)[i])^2 <=
#       4*sigma[i] + (L^{-1}u)[i]*((L^{-1}ytilde)[i]-(L^Tx)[i]);

subject to unlifted:
  sum{i in 1..n} x[i]*ytilde[i] - sum{j in 1..gdim} g[j]*y[j] <= 0;

#
# explicit forms of the complementarity constraints.
#  Note: it appears cplex doesn't handle this ampl syntax.
#
subject to comp_constraint_y{i in 1..gdim}:
  0 <= y[i] complements sum{j in 1..n} H[i,j]*x[j] >= g[i];

subject to comp_constraint_v{i in 1..n}:
  0 <= y[gdim+i] complements x[i] >= 0;

option ipopt_options "max_iter=300";
