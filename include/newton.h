# ifndef NEWTON_H
# define NEWTON_H

extern int NewtonMaxIter;  // default: 1000


// solves the equation func(x)=0  by Newton method
// places the root to xr
// acc is the maximal difference of func(xr) from 0
// starts iterations from xstart
// returns  0 if OK
//         -1 in case of zero derivative
//         -2 if iterations do not converge (derivative changes sign)
//         -3 if more than NewtonMaxIter iterations were made

int  NewtonSolve(float &xr,float acc,float func(float),float deriv(float),float xstart);


# endif
