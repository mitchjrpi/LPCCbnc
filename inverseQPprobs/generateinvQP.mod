#
# Follow construction of Jara-Moroni, Pang, Wachter: MP: DC for LPCCs
#
# Set up problems of the form:
# 
# min  || (x,b,c) - (xbar,bbar,cbar) ||_1
#  (x,b,c)
#
# s.t. -ux <= x <= ux, -ub <= b <= ub, -uc <= c <= uc
#       x in argmin{0.5 y^TQy + C^Ty : Ay >= b}
#              y
# 
# Equivalently:
# min         e^Tzx + e^Tzb + e^Tzc
#  (x,b,c,zx,zb,zc,lambda)
#
# s.t. -ux <= x <= ux, -ub <= b <= ub, -uc <= c <= uc
#        x+zx >=  xbar
#       -x+zx >= -xbar
#        b+zb >=  bbar
#       -b+zb >= -bbar
#        c+zc >=  cbar
#       -c+zc >= -cbar
#       -lambda >= -ulambda
#       Qx + c - A^T lambda >= 0
#      -Qx - c + A^T lambda >= 0
#      0 <= lambda \perp Ax - b >= 0
#
# A is m times n, everything else dimensioned appropriately.
# Q is positive definite.
#
# Q is a symmetric PD matrix with about 10 nonzeroes per row.
# Generate it by constructing a square matrix M and then Q=MM^T.
# M has nonzero diagonal entries plus two other nonzeroes per row,
# which should give about 10 nonzeroes per row of Q.
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
# Use the QPnonconvex/bilevel problems as an initial template.
#
# Generate everything as uniform01.
#  Include nonnegativity constraints in the inner problem.
#  Any w with Bw >= b and 0<=w<=u will then be feasible.
#
# Take e=vector of ones.
#
#param n > 0;   # dimension of x and c
#param m > 0;   # dimension of b and lambda

set Annz within {1..m,1..n} default {};
set Mnnz within {1..n,1..n} default {};
set Qnnz within {1..n,1..n} default {};

param A{Annz};
param M{Mnnz};
param Q{Qnnz};

param Arowsum{1..m};
param Arowcount{1..m};
param Acolcount{1..n};
param Qrowcount{1..n};

param xbb{j in 1..n} := Normal01();
param xbar{j in 1..n} := xbb[j] + Normal01();
param bbb{j in 1..m};
param bbar{j in 1..m};
param cbb{j in 1..n};
param cbar{j in 1..n};

param ux;
param uc;
param ub;
param ulambda;

param lambdabar {j in 1..m}:=Uniform(0,10);
param wbar {j in 1..m}:=Uniform(0,10);
param v{j in 1..m} := floor(Uniform(0,2)) binary;
param lambdabb {j in 1..m} := lambdabar[j]*v[j];
param wbb {j in 1..m} := wbar[j]*(1-v[j]);

#
param row;
param col;
param nz10 default 10;
param count;
param nnnz;
param nzcount;
#
var x{j in 1..n} >= -ux, <= ux;
var b{j in 1..m} >= -ub, <= ub;
var c{j in 1..n} >= -uc, <= uc;
var zx{j in 1..n};
var zb{j in 1..m};
var zc{j in 1..n};
var lambda{j in 1..m} >= 0, <= ulambda;
var wvar{j in 1..m};
#
var zcomp{1..m} binary;

#

minimize inverseQPobj: sum{i in 1..n} zx[i] + sum{j in 1..m} zb[j]
             + sum{k in 1..n} zc[k];

subject to gradcon1{i in 1..n}:
  c[i] + sum{j in 1..n : (i,j) in Qnnz} Q[i,j] * x[j]
       - sum{j in 1..m: (j,i) in Annz} A[j,i] * lambda[j] >= 0;

subject to gradcon2{i in 1..n}:
  - c[i] - sum{j in 1..n : (i,j) in Qnnz} Q[i,j] * x[j]
       + sum{j in 1..m: (j,i) in Annz} A[j,i] * lambda[j] >= 0;

subject to zxdef1{i in 1..n}:
  x[i] + zx[i] >= xbar[i];

subject to zxdef2{i in 1..n}:
  -x[i] + zx[i] >= -xbar[i];

subject to zbdef1{i in 1..m}:
  b[i] + zb[i] >= bbar[i];

subject to zbdef2{i in 1..m}:
  -b[i] + zb[i] >= -bbar[i];

subject to zcdef1{i in 1..n}:
  c[i] + zc[i] >= cbar[i];

subject to zcdef2{i in 1..n}:
  -c[i] + zc[i] >= -cbar[i];

subject to lowerfeas{i in 1..m}:
  sum{j in 1..n: (i,j) in Annz} A[i,j] * x[j] >= b[i];

# for MIP: enforce complementarity
subject to complambda{i in 1..m}:
  lambda[i] <= ulambda * zcomp[i];

subject to compw{i in 1..m}:
  sum{j in 1..n : (i,j) in Annz} A[i,j] * x[j] - b[i]
    <= (ux * Arowsum[i] + ub) * (1-zcomp[i]);

subject to calcw{i in 1..m}:
  wvar[i] = sum{j in 1..n : (i,j) in Annz} A[i,j] * x[j] - b[i];

#
# explicit forms of the complementarity constraints.
#  Note: it appears cplex doesn't handle this ampl syntax.
#
# subject to comp_constraint_y{i in 1..gdim}:
#   0 <= y[i] complements sum{j in 1..n} H[i,j]*x[j] >= g[i];
# 
# subject to comp_constraint_v{i in 1..n}:
#   0 <= y[gdim+i] complements x[i] >= 0;
