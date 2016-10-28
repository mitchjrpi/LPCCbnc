Solving problems of the form:

min  c^Tx  +  d^Ty
st    Ax   +  By  >=  b
     0 <= y \perp w := q + Nx + My >= 0

c is an n-vector.
d is an m-vector.
f is a  k-vector.
All other vectors and matrices dimensioned appropriately.

format_3 is a sparse format.

Row 1:   [n,m,k]
Row 2:   [c(1),...,c(n)]
Row 3:   [d(1),...,d(m)]
Row 4:   [b(1),...,b(k)]
Row 5:   [q(1),...,q(m)]
Rows 6 to 11: A as:
  row 6:  [[k,n,#nonzeroes in A],
  row 7:   [0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
                 #nonzeroes in rows 1 to 3,...,#nonzeroes in rows 1 to k-1],
  row 8:   [#nonzeroes in row 1,#nonzeroes in row 2,...,#nonzeroes in row k],
  row 9:   [column indices of the nonzeroes in A],
  row 10:  [values of the nonzeroes in A]
  row 11: ]
Rows 12 to 17: B as:
  row 12: [[k,m,#nonzeroes in B],
  row 13:  [0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
                #nonzeroes in rows 1 to 3,...,#nonzeroes in rows 1 to k-1],
  row 14:  [#nonzeroes in row 1,#nonzeroes in row 2,...,#nonzeroes in row k],
  row 15:  [column indices of the nonzeroes in B],
  row 16:  [values of the nonzeroes in B]
  row 17: ]
Rows 18 to 23: N as:
  row 18: [[m,n,#nonzeroes in N],
  row 19:  [0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
                #nonzeroes in rows 1 to 3,...,#nonzeroes in rows 1 to m-1],
  row 20:  [#nonzeroes in row 1,#nonzeroes in row 2,...,#nonzeroes in row m],
  row 21:  [column indices of the nonzeroes in N],
  row 22:  [values of the nonzeroes in N]
  row 23: ]
Rows 24 to 29: M as:
  row 24: [[m,m,#nonzeroes in M],
  row 25:  [0,#nonzeroes in row 1,#nonzeroes in rows 1 and 2,
                #nonzeroes in rows 1 to 3,...,#nonzeroes in rows 1 to m-1],
  row 26:  [#nonzeroes in row 1,#nonzeroes in row 2,...,#nonzeroes in row m],
  row 27:  [column indices of the nonzeroes in M],
  row 28:  [values of the nonzeroes in M]
  row 29: ]
