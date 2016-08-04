Set up problems of the form:

min  c^Tx + d^Tw
st   Ax + Bw >= b
     0 <= w <= u
     x in argmin_x { 0.5 x^T Qx + w^Tx : Hx >= g}.

Equivalently:
min  c^Tx + d^Tw
st   Ax + Bw >= b
     0 <= w <= u
     Qx + w - H^Ty = 0
     0 <= y \perp Hx - g >= 0.

Construct the instances in format_3.

The matrices are stored in a sparse manner:
[[#rows,#cols,#nonzeroes],
[0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
#nonzeroes in rows 1 and 2 and 3,...,
#nonzeroes in rows 1..(#rows-1)],
[#nonzeroes in row 1,#nonzeroes in row 2,
#nonzeroes in row 3,...,
#nonzeroes in row (#rows-1),#nonzeroes in last row],
[(column of each nonzero) - 1],
[the nonzeroes]
]

For our problem format, the order of storage of the data
is as follows, in 29 lines:
 
line 1:
[dim x + dim w, dim y, dim b + 3*dim w (for bounds on w and the
gradient KKT condition for the inner problem)]

line 2:
[coeffs of c and d]

line 3:
[zeroes for the coeffs of y in objective function]

line 4:
[coeffs of b, zeroes, -components of u, zeroes]

line 5:
[-components of g]

lines 6-11:
[(store the A,B, the I and -I for the bounds on w,
and also the matrices Q and I for Qx+w-H^Ty=0)
[dim b + 3*dim w,dim x + dim w,#nonzeroes in A and B + 2*dim w],
[dim b + 3*dim w entries giving the cumulative sum of nonzeroes
in the rows, with an extra 0 at the beginning, and no entry for
the final row],
[dim b + 3*dim w entries giving the numbers of nonzeroes in each row],
[(column entry of each nonzero) - 1],
[actual nonzero values]
]

lines 12-17:
[(store the matrix of zeroes for the coefficients of y in the outer constraints
and also H^T from the gradient condition)
[dim b + 3*dim w,dim y,#nonzeroes in H^T],
[dim b + 2*dim w + 1 zeroes, cumulative sums of nonzeroes in rows
of H^T (skipping the final row)],
[dim b + 2*dim w zeroes, #nonzeroes in each row of H^T],
[(column of each nonzero in H^T) - 1],
[nonzeroes in H^T]
]

lines 18-23:
[(store H and a matrix of zeroes for the coefficients of w in the
complementarity constraints)
[dim y,dim x + dim w,#nonzeroes in H],
[0, cumulative sums of nonzeroes in rows of H (skipping the final row)],
[#nonzeroes in each row of H],
[(column of each nonzero in H) - 1],
[nonzeroes in H]
]

lines 24-29:
[(store a matrix of zeroes for the coefficients of y on the
right hand side of the complementarity constraints)
[dim y,dim y,0],
[dim y zeroes],
[dim y zeroes],
[],
[]
]
