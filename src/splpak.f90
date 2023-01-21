!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!

    module splpak_module

   ! implicit none

    contains

subroutine bascmpd(x,nderiv,xmin,nodes,icol,basm)
double precision x
double precision xmin
double precision basm
double precision dx
double precision dxin
double precision xb
double precision bas1
double precision z
double precision fact
double precision z1
!
!  This routine does basis function computations for natural
!  splines.  This routine is called by routines SPLCW and SPLDE
!  to compute ICOL and BASM, which are defined as follows:
!
!     The MDIM indices in IB (defined through common) determine
!     a specific node in the node grid (see routine SPLCC for a
!     description of the node grid).  Every node is associated
!     with an MDIM-dimensional basis function and a corresponding
!     column in the least squares matrix (or element of the
!     coefficient vector).  The column index (which may be thought
!     of as a linear address for the MDIM-dimensional node grid)
!     corresponding to the specified node is computed as ICOL.  The
!     associated basis function evaluated at X (an MDIM-vector) is
!     computed as BASM (a scalar).
!
!  In case NDERIV is not all zero, BASM will not be the value of
!  the basis function but rather a partial derivative of that
!  function as follows:
!
!     The order of the partial derivative in the direction of the
!     IDIM coordinate is NDERIV(IDIM) (for IDIM <= MDIM).  This
!     routine will compute incorrect values if NDERIV(IDIM) is not
!     in the range 0 to 2.
!
!
dimension x(4),nderiv(4),xmin(4),nodes(4)
!
!  The technique of this routine is to transform the independent
!  variable in each dimension such that the nodes fall on
!  suitably chosen integers.  On this transformed space, the
!  1-dimensional basis functions and their derivatives have a
!  particularly simple form.  The desired MDIM-dimensional basis
!  function (or any of its partial derivatives) is computed as
!  a product of such 1-dimensional functions (tensor product
!  method of defining multi-dimensional splines).  The values
!  which determine the location of the nodes, and hence the
!  above transform, are passed through common and the argument
!  list.
!
common /splcomd/dx(4),dxin(4),mdim,ib(4),ibmn(4),ibmx(4)
save
!
!  ICOL will be a linear address corresponding to the indices in IB.
!
icol = 0
!
!  BASM will be M-dimensional basis function evaluated at X.
!
basm = 1.d0
do 121 idim = 1,mdim
!
!  Compute ICOL by Horner's method.
!
    mdmid = mdim + 1 - idim
    icol = nodes(mdmid)*icol + ib(mdmid)
!
!  NGO depends upon function type and NDERIV.
!
    ntyp = 1
!
!  Function type 1 (left linear) for IB = 0 or 1.
!
    if (ib(idim)<=1) go to 101
    ntyp = 2
!
!  Function type 2 (chapeau function) for 2 LT IB LT NODES-2.
!
    if (ib(idim)<nodes(idim)-2) go to 101
    ntyp = 3
!
!  Function type 3 (right linear) for IB = NODES-2 or NODES-1.
!
101     ngo = 3*ntyp + nderiv(idim) - 2
!
!  XB is X value of node IB (center of basis function).
!
    xb = xmin(idim) + dble(ib(idim))*dx(idim)
!
!  BAS1 will be the 1-dimensional basis function evaluated at X.
!
    bas1 = 0.d0
    go to (102,103,104,105,106,108,110,113,117) ngo
!
!  Function type 1 (left linear) is mirror image of function type 3.
!
!  Transform so that XB is at 2 and the other nodes are at the integers
!  (with ordering reversed to form a mirror image).
!
102     z = dxin(idim)* (xb-x(idim)) + 2.d0
    go to 111
!
!  1st derivative.
!
103     fact = -dxin(idim)
    go to 114
!
!  2nd derivative.
!
104     fact = -dxin(idim)
    go to 118
!
!  Function type 2 (chapeau function).
!
!  Transform so that XB is at the origin and the other nodes are at
!  the integers.
!
105     z = abs(dxin(idim)* (x(idim)-xb)) - 2.d0
!
!  This chapeau function is then that unique cubic spline which is
!  identically zero for ABS(Z) GE 2 and is 1 at the origin.  This
!  function is the general interior node basis function.
!
    if (z>=0.d0) go to 120
    bas1 = -.25d0*z**3
    z = z + 1.d0
    if (z>=0.d0) go to 120
    bas1 = bas1 + z**3
    go to 120
!
!  1st derivative.
!
106     z = x(idim) - xb
    fact = dxin(idim)
    if (z<0.d0) fact = -fact
    z = fact*z - 2.d0
    if (z>=0.d0) go to 120
    bas1 = -.75d0*z**2
    z = z + 1.d0
    if (z>=0.d0) go to 107
    bas1 = bas1 + 3.d0*z**2
107     bas1 = fact*bas1
    go to 120
!
!  2nd derivative.
!
108     fact = dxin(idim)
    z = fact*abs(x(idim)-xb) - 2.d0
    if (z>=0.d0) go to 120
    bas1 = -1.5d0*z
    z = z + 1.d0
    if (z>=0.d0) go to 109
    bas1 = bas1 + 6.d0*z
109     bas1 = (fact**2)*bas1
    go to 120
!
!  Function type 3 (right linear).
!
!  Transform so that XB is at 2 and the other nodes are at the integers.
!
110     z = dxin(idim)* (x(idim)-xb) + 2.d0
!
!  This right linear function is defined to be that unique cubic spline
!  which is identically zero for Z <= 0 and is 3*Z-3 for Z GE 2.
!  This function (obviously having zero 2nd derivative for
!  Z GE 2) is used for the two nodes nearest an edge in order
!  to generate natural splines, which must by definition have
!  zero 2nd derivative at the boundary.
!
!  Note that this method of generating natural splines also provides
!  a linear extrapolation which has 2nd order continuity with
!  the interior splines at the boundary.
!
111     if (z<=0.d0) go to 120
    if (z>=2.d0) go to 112
    bas1 = .5d0*z**3
    z = z - 1.d0
    if (z<=0.d0) go to 120
    bas1 = bas1 - z**3
    go to 120
112     bas1 = 3.d0*z - 3.d0
    go to 120
!
!  1st derivative.
!
113     fact = dxin(idim)
114     z = fact* (x(idim)-xb) + 2.d0
    if (z<=0.d0) go to 120
    if (z>=2.d0) go to 116
    bas1 = 1.5d0*z**2
    z = z - 1.d0
    if (z<=0.d0) go to 115
    bas1 = bas1 - 3.d0*z**2
115     bas1 = fact*bas1
    go to 120
116     bas1 = 3.d0*fact
    go to 120
!
!  2nd derivative.
!
117     fact = dxin(idim)
118     z = fact* (x(idim)-xb) + 2.d0
    z1 = z - 1.d0
    if (abs(z1)>=1.d0) go to 120
    bas1 = 3.d0*z
    if (z1<=0.d0) go to 119
    bas1 = bas1 - 6.d0*z1
119     bas1 = (fact**2)*bas1
120     basm = basm*bas1
121 continue
icol = icol + 1
return
end
!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!
! SUBROUTINE CFAERR (IERR,MESS,LMESS)
!
! PURPOSE        To print an error number and an error message
!                or just an error message.
!
! USAGE          CALL CFAERR (IERR,MESS,LMESS)
!
! ARGUMENTS
! ON INPUT       IERR
!                  The error number (printed only if non-zero).
!
!                MESS
!                  Message to be printed.
!
!                LMESS
!                  Number of characters in mess (<= 130).
!
! ARGUMENTS
! ON OUTPUT      None
!
! I/O            The message is writen to unit 6.
!
! ******************************************************************
!
subroutine cfaerr (ierr,mess,lmess)
!
character *(*) mess
!
if (ierr /= 0) write (6,'(A,I5)') ' IERR=', ierr
write (6,'(A)') mess(1:lmess)
!
return
end
! PACKAGE SPLPAK         Documentation for user entries follows
!                        the general package information.
!
! LATEST REVISION        August, 1998
!
! PURPOSE                This package contains routines for fitting
!                        (least squares) a multidimensional cubic spline
!                        to arbitrarily located data.  It also contains
!                        routines for evaluating this spline (or its
!                        partial derivatives) at any point.
!
!                        Coefficient calculation is performed in
!                        subroutines SPLCC or SPLCW and evaluation is
!                        performed by functions SPLFE or SPLDE.
!
! USAGE                  Package SPLPAK contains four user entries --
!                        SPLCC, SPLCW, SPLFE, AND SPLDE.
!
!                        The user first calls SPLCC by
!
!                          CALL SPLCC (NDIM,XDATA,L1XDAT,YDATA,NDATA,
!                                      XMIN,XMAX,NODES,XTRAP,COEF,NCF,
!                                      WORK,NWRK,IERROR)
!
!                        or SPLCW by
!
!                          CALL SPLCW (NDIM,XDATA,L1XDATA,YDATA,WDATA,
!                                      NDATA,XMIN,XMAX,NODES,XTRAP,
!                                      COEF,NCF,WORK,NWRK,IERROR)
!
!                        The parameter NDATA in the call to SPLCW
!                        enables the user to weight some of the data
!                        points more heavily than others.  Both
!                        routines return a set of coefficients in the
!                        array COEF.  These coefficients are
!                        subsequently used in the computation of
!                        function values and partial derivatives.
!                        To compute values on the spline approximation
!                        the user then calls SPLFE or SPLDE any
!                        number of times in any order provided that
!                        the values of the inputs, NDIM, COEF, XMIN,
!                        XMAX, and NODES, are preserved between calls.
!
!*PL*ERROR* Comment line too long
!                        SPLFE and SPLDE are called in the following way:
!
!                          F = SPLFE (NDIM,X,COEF,XMIN,XMAX,NODES,
!                                     IERROR)
!
!                        or
!
!                          F = SPLDE (NDIM,X,NDERIV,COEF,XMIN,XMAX,
!                                     NODES,IERROR)
!
!                        The routine SPLFE returns an interpolated
!                        value at the point defined by the array X.
!                        SPLDE affords the user the additional
!                        capability of calculating an interpolated
!                        value for one of several partial derivatives
!                        specified by the array NDERIV.
!
! I/O                    None, except for error messages printed by
!                        calls to CFSARR, if an error is detected.
!
! PRECISION              Single
!
! REQUIRED LIBRARY       SUPRLS, CFAERR
! FILES
!
! LANGUAGE               FORTRAN
!
! HISTORY                Developed in 1972-73 by NCAR's
!                        Scientific Computing Division.
!
!                        Cleaned up and added to the Ngmath library in
!                        1998.
!
! PORTABILITY            FORTRAN 77
!
!***********************************************************************
!
! SUBROUTINE SPLCCD(NDIM,XDATA,L1XDAT,YDATA,NDATA,XMIN,XMAX,NODES,
!                   XTRAP,COEF,NCF,WORK,NWRK,IERROR)
!
! DIMENSION OF           XDATA(NDATA,L1XDAT),YDATA(NDATA),XMIN(NDIM),
! ARGUMENTS              XMAX(NDIM),NODES(NDIM),COEF(NCF),WORK(NWRK)
!
! PURPOSE                N-dimensional cubic spline coefficient
!                        calculation by least squares.
!
! USAGE                  The usage and arguments of this routine are
!                        identical to those for SPLCW except for the
!                        omission of the array of weights, WDATA.  See
!                        entry SPLCW description immediately below for a
!                        complete description.
!
!                        CALL SPLCCD(NDIM,XDATA,L1XDAT,YDATA,NDATA,XMIN,
!                                    XMAX,NODES,XTRAP,COEF,NCF,WORK,
!                                    NWRK,IERROR)
!
!***********************************************************************
!
! SUBROUTINE SPLCWD(NDIM,XDATA,L1XDAT,YDATA,WDATA,NDATA,XMIN,XMAX,
!                   NODES,XTRAP,COEF,NCF,WORK,NWRK,IERROR)
!
!
! DIMENSION OF           XDATA(L1XDAT,NDATA),YDATA(NDATA),WDATA(NDATA),
! ARGUMENTS              XMIN(NDIM),XMAX(NDIM),NODES(NDIM),COEF(NCF),
!                        WORK(NWRK)
!
! PURPOSE                N-dimensional cubic spline coefficient
!                        calculation by weighted least squares on
!                        arbitrarily located data.
!
!                        A grid of evenly spaced nodes in NDIM space is
!                        defined by the arguments XMIN, XMAX and NODES.
!                        A linear basis for the class of natural splines
!                        on these nodes is formed, and a set of
!                        corresponding coefficients is computed in the
!                        array COEF.  These coefficients are chosen to
!                        minimize the weighted sum of squared errors
!                        between the spline and the arbitrarily located
!                        data values described by the arguments XDATA,
!                        YDATA and NDATA.  The smoothness of the spline
!                        in data sparse areas is controlled by the
!                        argument XTRAP.
!
! NOTE                   In order to understand the arguments of this
!                        routine, one should realize that the node grid
!                        need not bear any particular relation to the
!                        data points.  In the theory of exact-fit
!                        interpolatory splines, the nodes would in fact
!                        be data locations, but in this case they serve
!                        only to define the class of splines from which
!                        the approximating function is chosen.  This
!                        node grid is a rectangular arrangement of
!                        points in NDIM space, with the restriction that
!                        along any coordinate direction the nodes are
!                        equally spaced.  The class of natural splines
!                        on this grid of nodes (NDIM-cubic splines whose
!                        2nd derivatives normal to the boundaries are 0)
!                        has as many degrees of freedom as the grid has
!                        nodes.  Thus the smoothness or flexibility of
!                        the splines is determined by the choice of the
!                        node grid.
!
! USAGE                  CALL SPLCW (NDIM,XDATA,L1XDAT,YDATA,WDATA,
!                                    NDATA,XMIN,XMAX,NODES,XTRAP,COEF,
!                                    NCF,WORK,NWRK,IERROR)
!
!                        The spline (or its derivatives) may then be
!                        evaluated by using function SPLFE (or SPLDE).
!
! ARGUMENTS
!
! ON INPUT               NDIM
!                          The dimensionality of the problem.  The
!                          spline is a function of NDIM variables or
!                          coordinates and thus a point in the
!                          independent variable space is an NDIM vector.
!                          NDIM must be in the range 1 <= NDIM <= 4.
!
!                        XDATA
!                          A collection of locations for the data
!                          values, i.e., points from the independent
!                          variable space.  This collection is a
!                          2-dimensional array whose 1st dimension
!                          indexes the NDIM coordinates of a given point
!                          and whose 2nd dimension labels the data
!                          point.  For example, the data point with
!                          label IDATA is located at the point
!                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) where
!                          the elements of this vector are the values of
!                          the NDIM coordinates.  The location, number
!                          and ordering of the data points is arbitrary.
!                          The dimension of XDATA is assumed to be
!                          XDATA(L1XDAT,NDATA).
!
!                        L1XDAT
!                          The length of the 1st dimension of XDATA in
!                          the calling program.  L1XDAT must be >=
!                          NDIM.
!
!                               NOTE:  For 1-dimensional problems L1XDAT
!                                      is usually 1.
!
!                        YDATA
!                          A collection of data values corresponding to
!                          the points in XDATA.  YDATA(IDATA) is the
!                          data value associated with the point
!                          (XDATA(1,IDATA),...,XDATA(NDIM,IDATA)) in the
!                          independent variable space.  The spline whose
!                          coefficients are computed by this routine
!                          approximates these data values in the least
!*PL*ERROR* Comment line too long
!                          squares sense.  The dimension is assumed to be
!                          YDATA(NDATA).
!
!                        WDATA
!                          A collection of weights.  WDATA(IDATA) is a
!                          weight associated with the data point
!                          labelled IDATA.  It should be non-negative,
!                          but may be of any magnitude.  The weights
!                          have the effect of forcing greater or lesser
!                          accuracy at a given point as follows: this
!                          routine chooses coefficients to minimize the
!                          sum over all data points of the quantity
!
!                            (WDATA(IDATA)*(YDATA(IDATA) - spline value
!                            at XDATA(IDATA)))**2.
!
!                          Thus, if the reliability
!                          of a data point is known to be low, the
!                          corresponding weight may be made small
!                          (relative to the other weights) so that the
!                          sum over all data points is affected less by
!                          discrepencies at the unreliable point.  Data
!                          points with zero weight are completely
!                          ignored.
!
!                               NOTE:  If WDATA(1) is < 0, the other
!                                      elements of WDATA are not
!                                      referenced, and all weights are
!                                      assumed to be unity.
!
!                          The dimension is assumed to be WDATA(NDATA)
!                          unless WDATA(1) < 0., in which case the
!                          dimension is assumed to be 1.
!
!                        NDATA
!                          The number of data points mentioned in the
!                          above arguments.
!
!                        XMIN
!                          A vector describing the lower extreme corner
!                          of the node grid.  A set of evenly spaced
!                          nodes is formed along each coordinate axis
!                          and XMIN(IDIM) is the location of the first
!                          node along the IDIM axis.  The dimension is
!                          assumed to be XMIN(NDIM).
!
!                        XMAX
!                          A vector describing the upper extreme corner
!                          of the node grid.  A set of evenly spaced
!                          nodes is formed along each coordinate axis
!                          and XMAX(IDIM) is the location of the last
!                          node along the IDIM axis.  The dimension is
!                          assumed to be XMAX(NDIM).
!
!                        NODES
!                          A vector of integers describing the number of
!                          nodes along each axis.  NODES(IDIM) is the
!                          number of nodes (counting endpoints) along
!                          the IDIM axis and determines the flexibility
!                          of the spline in that coordinate direction.
!                          NODES(IDIM) must be >= 4, but may be as
!                          large as the arrays COEF and WORK allow.
!                          The dimension is assumed to be NODES(NDIM).
!
!                          NOTE:  The node grid is completely defined by
!                                 the arguments XMIN, XMAX and NODES.
!                                 The spacing of this grid in the IDIM
!                                 coordinate direction is:
!
!                                   DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) /
!                                              (NODES(IDIM)-1).
!
!                                 A node in this grid may be indexed by
!                                 an NDIM vector of integers
!                                 (IN(1),...,IN(NDIM)) where
!                                 1 <= IN(IDIM) <= NODES(IDIM).
!                                 The location of such a node may be
!                                 represented by an NDIM vector
!                                 (X(1),...,X(NDIM)) where
!                                 X(IDIM) = XMIN(IDIM) + (IN(IDIM)-1) *
!                                                         DX(IDIM).
!
!                        XTRAP
!                          A parameter to control extrapolation to data
!                          sparse areas.  The region described by XMIN
!                          and XMAX is divided into rectangles, the
!                          number of which is determined by NODES, and
!                          any rectangle containing a disproportionately
!                          small number of data points is considered to
!                          be data sparse (rectangle is used here to
!                          mean NDIM-dimensional rectangle).  If XTRAP
!                          is nonzero the least squares problem is
!                          augmented with derivative constraints in the
!                          data sparse areas to prevent the matrix from
!                          becoming poorly conditioned.  XTRAP serves as
!                          a weight for these constraints, and thus may
!                          be used to control smoothness in data sparse
!                          areas.  Experience indicates that unity is a
!                          good first guess for this parameter.
!
!                               NOTE:  If XTRAP is zero, substantial
!                                      portions of the routine will be
!                                      skipped, but a singular matrix
!                                      can result if large portions of
!                                      the region are without data.
!
!                        NCF
!                          The length of the array COEF in the calling
!                          program.  If NCF is <
!                          NODES(1)*...*NODES(NDIM), a fatal error is
!                          diagnosed.
!
!                        WORK
!                          A workspace array for solving the least
!                          squares matrix generated by this routine.
!                          Its required size is a function of the total
!                          number of nodes in the node grid.  This
!                          total, NCOL = NODES(1)*...*NODES(NDIM), is
!                          also the number of columns in the least
!                          squares matrix.  The length of the array WORK
!                          must equal or exceed NCOL*(NCOL+1).
!
!                        NWRK
!                          The length of the array WORK in the calling
!                          program.  If
!                          NCOL = NODES(1)*...*NODES(NDIM) is the total
!                          number of nodes, then a fatal error is
!                          diagnosed if NWRK is less than
!                          NCOL*(NCOL+1).
!
! ON OUTPUT              COEF
!                          The array of coefficients computed by this
!                          routine.  Each coefficient corresponds to a
!                          particular basis function which in turn
!                          corresponds to a node in the node grid.  This
!                          correspondence between the node grid and the
!                          array COEF is as if COEF were an
!                          NDIM-dimensional Fortran array with
!                          dimensions NODES(1),...,NODES(NDIM), i.e., to
!                          store the array linearly, the leftmost
!                          indices are incremented most frequently.
!                          Hence the length of the COEF array must equal
!                          or exceed the total number of nodes, which is
!                          NODES(1)*...*NODES(NDIM).  The computed array
!                          COEF may be used with function SPLFE
!                          (or SPLDE) to evaluate the spline (or its
!                          derivatives) at an arbitrary point in NDIM
!*PL*ERROR* Comment line too long
!                          space.  The dimension is assumed to be COEF(NCF).
!
!                        WORK
!                          The workspace containing intermediate
!                          calculations.  It need not be saved.
!
!                        IERROR
!                          An error flag with the following meanings:
!                              0  No error.
!                            101  NDIM is < 1 or is > 4.
!                            102  NODES(IDIM) is < 4 fOR some IDIM.
!                            103  XMIN(IDIM) = XMAX(IDIM) for some IDIM.
!                            104  NCF (size of COEF) is
!                                 < NODES(1)*...*NODES(NDIM).
!                            105  NDATA is < 1.
!                            106  NWRK (size of WORK) is too small.
!                            107  SUPRLS failure (usually insufficient
!                                 data) -- ordinarily occurs only if
!                                 XTRAP is zero or WDATA contains all
!                                 zeros.
!
! ALGORITHM              An overdetermined system of linear equations
!                        is formed -- one equation for each data point
!                        plus equations for derivative constraints.
!                        This system is solved using subroutine SUPRLSD.
!
! ACCURACY               If there is exactly one data point in the
!                        near vicinity of each node and no extra data,
!                        the resulting spline will agree with the
!                        data values to machine accuracy.  However, if
!                        the problem is overdetermined or the sparse
!                        data option is utilized, the accuracy is hard
!                        to predict.  Basically, smooth functions
!                        require fewer nodes than rough ones for the
!                        same accuracy.
!
! TIMING                 The execution time is roughly proportional
!                        to NDATA*NCOF**2 where NCOF = NODES(1)*...*
!                        NODES(NDIM).
!
!***********************************************************************
!
subroutine splccd(ndim,xdata,l1xdat,ydata,ndata,xmin,xmax,nodes, &
                  xtrap,coef,ncf,work,nwrk,ierror)
double precision xdata
double precision ydata
double precision xmin
double precision xmax
double precision xtrap
double precision coef
double precision work
double precision w
dimension xdata(l1xdat,ndata),ydata(ndata),xmin(ndim),xmax(ndim), &
          nodes(ndim),coef(ncf),work(nwrk)
dimension w(1)
save
!
w(1) = -1.d0
call splcwd(ndim,xdata,l1xdat,ydata,w,ndata,xmin,xmax,nodes,xtrap, &
            coef,ncf,work,nwrk,ierror)
!
return
end
!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!
subroutine splcwd(ndim,xdata,l1xdat,ydata,wdata,ndata,xmin,xmax, &
                 nodes,xtrap,coef,ncf,work,nwrk,ierror)
double precision xdata
double precision ydata
double precision wdata
double precision xmin
double precision xmax
double precision xtrap
double precision coef
double precision work
double precision x
double precision dx
double precision dxin
double precision spcrit
double precision xrng
double precision swght
double precision rowwt
double precision rhs
double precision basm
double precision reserr
double precision totlwt
double precision bump
double precision wtprrc
double precision expect
double precision dcwght
dimension xdata(l1xdat,ndata),ydata(ndata),wdata(ndata), &
          xmin(ndim),xmax(ndim),nodes(ndim),coef(ncf),work(nwrk)
dimension x(4),nderiv(4),in(4),inmx(4)
common /splcomd/dx(4),dxin(4),mdim,ib(4),ibmn(4),ibmx(4)
save
!
!  The restriction that NDIM be less than are equat to 4 can be
!  eliminated by increasing the above dimensions, but the required
!  length of WORK becomes quite large.
!
!  SPCRIT is used to determine data sparseness as follows -
!  the weights assigned to all data points are totaled into the
!  variable TOTLWT. (If no weights are entered, it is set to
!  NDATA.)  Each node of the node network is assigned a
!  rectangle (in which it is contained) and the weights of all
!  data points which fall in that rectangle are totaled.  If that
!  total is less than SPCRIT*EXPECT (EXPECT is defined below),
!  then the node is ascertained to be in a data sparse location.
!  EXPECT is that fraction of TOTLWT that would be expected by
!  comparing the area of the rectangle with the total area under
!  consideration.
!
data spcrit/.75d0/
!
ierror = 0
mdim = ndim
if (mdim<1 .or. mdim>4) go to 127
ncol = 1
do 101 idim = 1,mdim
    nod = nodes(idim)
    if (nod<4) go to 128
!
!  Number of columns in least squares matrix = number of coefficients =
!  product of nodes over all dimensions.
!
    ncol = ncol*nod
    xrng = xmax(idim) - xmin(idim)
    if (xrng==0.d0) go to 129
!
!  DX(IDIM) is the node spacing along the IDIM coordinate.
!
    dx(idim) = xrng/dble(nod-1)
    dxin(idim) = 1.d0/dx(idim)
    nderiv(idim) = 0
101 continue
if (ncol>ncf) go to 130
nwrk1 = 1
mdata = ndata
if (mdata<1) go to 131
!
!  SWGHT is a local variable = XTRAP, and can be considered a smoothing
!  weight for data sparse areas.  If SWGHT == 0, no smoothing
!  computations are performed.
!
swght = xtrap
!
!  Set aside workspace for counting data points.
!
if (swght/=0.d0) nwrk1 = ncol + 1
!
!  NWLFT is the length of the remaining workspace.
!
nwlft = nwrk - nwrk1 + 1
if (nwlft<1) go to 132
irow = 0
!
!  ROWWT is used to weight rows of the least squares matrix.
!
rowwt = 1.d0
!
!  Loop through all data points, computing a row for each.
!
do 108 idata = 1,mdata
!
!  WDATA(1)<0 means weights have not been entered.  In that case,
!  ROWWT is left equal to  1. for all points.  Otherwise ROWWT is
!  equal to WDATA(IDATA).
!
!  Every element of the row, as well as the corresponding right hand
!  side, is multiplied by ROWWT.
!
    if (wdata(1)<0.d0) go to 102
    rowwt = wdata(idata)
!
!  Data points with 0 weight are ignored.
!
    if (rowwt==0.d0) go to 108
102     irow = irow + 1
!
!  One row of the least squares matrix corresponds to each data
!  point.  The right hand for that row will correspond to the
!  function value YDATA at that point.
!
    rhs = rowwt*ydata(idata)
    do 103 idim = 1,mdim
        x(idim) = xdata(idim,idata)
103     continue
!
!  The COEF array serves as a row of least squares matrix.
!  Its value is zero except for columns corresponding to functions
!  which are nonzero at X.
!
    do 104 icol = 1,ncol
        coef(icol) = 0.d0
104     continue
!
!  Compute the indices of basis functions which are nonzero at X.
!  IBMN is in the range 0 to nodes-2 and IBMX is in range 1
!  to NODES-1.
!
    do 105 idim = 1,mdim
        nod = nodes(idim)
        it = dxin(idim)* (x(idim)-xmin(idim))
        ibmn(idim) = min0(max0(it-1,0),nod-2)
        ib(idim) = ibmn(idim)
        ibmx(idim) = max0(min0(it+2,nod-1),1)
105     continue
!
!  Begining of basis index loop - traverse all indices corresponding
!  to basis functions which are nonzero at X.  The indices are in
!  IB and are passed through common to BASCMP.
!
106     call bascmpd(x,nderiv,xmin,nodes,icol,basm)
!
!  BASCMP computes ICOL and BASM where BASM is the value at X of
!  the N-dimensional basis function corresponding to column ICOL.
!
    coef(icol) = rowwt*basm
!
!  Increment the basis indices.
!
    do 107 idim = 1,mdim
        ib(idim) = ib(idim) + 1
        if (ib(idim)<=ibmx(idim)) go to 106
        ib(idim) = ibmn(idim)
107     continue
!
!  End of basis index loop.
!
!
!  Send a row of the least squares matrix to the reduction routine.
!
    call suprld(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef,reserr, &
                lserr)
    if (lserr/=0) go to 133
108 continue
!
!  Row computations for all data points are now complete.
!
!  If SWGHT==0, the least squares matrix is complete and no
!  smoothing rows are computed.
!
if (swght==0.d0) go to 126
!
!  Initialize smoothing computations for data sparse areas.
!  Derivative constraints will always have zero right hand side.
!
rhs = 0.d0
nrect = 1
!
!  Initialize the node indices and compute number of rectangles
!  formed by the node network.
!
do 109 idim = 1,mdim
    in(idim) = 0
    inmx(idim) = nodes(idim) - 1
    nrect = nrect*inmx(idim)
109 continue
!
!  Every node is assigned an element of the workspace (set aside
!  previously) in which data points are counted.
!
do 110 iin = 1,ncol
    work(iin) = 0.d0
110 continue
!
!  Assign each data point to a node, total the assignments for
!  each node, and save in the workspace.
!
totlwt = 0.d0
do 112 idata = 1,mdata
!
!  BUMP is the weight associated with the data point.
!
    bump = 1.d0
    if (wdata(1)>=0.d0) bump = wdata(idata)
    if (bump==0.d0) go to 112
!
!  Find the nearest node.
!
    iin = 0
    do 111 idimc = 1,mdim
        idim = mdim + 1 - idimc
        inidim = int(dxin(idim)* (xdata(idim,idata)-xmin(idim))+ &
                 .5d0)
!
!  Points not in range (+ or - 1/2 node spacing) are not counted.
!
        if (inidim<0 .or. inidim>inmx(idim)) go to 112
!
!  Compute linear address of node in workspace by Horner's method.
!
        iin = (inmx(idim)+1)*iin + inidim
111     continue
!
!  Bump counter for that node.
!
    work(iin+1) = work(iin+1) + bump
    totlwt = totlwt + bump
112 continue
!
!  Compute the expected weight per rectangle.
!
wtprrc = totlwt/dble(nrect)
!
!  IN contains indices of the node (previously initialized).
!  IIN will be the linear address of the node in the workspace.
!
iin = 0
!
!  Loop through all nodes, computing derivative constraint rows
!  for those in data sparse locations.
!
!  Begining of node index loop - traverse all node indices.
!  The indices are in IN.
!
113 iin = iin + 1
expect = wtprrc
!
!  Rectangles at edge of network are smaller and hence less weight
!  should be expected.
!
do 114 idim = 1,mdim
    if (in(idim)==0 .or. in(idim)==inmx(idim)) expect = .5d0* &
        expect
114 continue
!
!  The expected weight minus the actual weight serves to define
!  data sparseness and is also used to weight the derivative
!  constraint rows.
!
!  There is no constraint if not data sparse.
!
if (work(iin)>=spcrit*expect) go to 124
dcwght = expect - work(iin)
do 115 idim = 1,mdim
    inidim = in(idim)
!
!  Compute the location of the node.
!
    x(idim) = xmin(idim) + dble(inidim)*dx(idim)
!
!  Compute the indices of the basis functions which are non-zero
!  at the node.
!
    ibmn(idim) = inidim - 1
    ibmx(idim) = inidim + 1
!
!  Distinguish the boundaries.
!
    if (inidim==0) ibmn(idim) = 0
    if (inidim==inmx(idim)) ibmx(idim) = inmx(idim)
!
!  Initialize the basis indices.
!
    ib(idim) = ibmn(idim)
115 continue
!
!  Multiply by the extrapolation parameter (this acts as a
!  smoothing weight).
!
dcwght = swght*dcwght
!
!  The COEF array serves as a row of the least squares matrix.
!  Its value is zero except for columns corresponding to functions
!  which are non-zero at the node.
!
do 116 icol = 1,ncol
    coef(icol) = 0.d0
116 continue
!
!  The 2nd derivative of a function of MDIM variables may be thought
!  of as a symmetric MDIM x MDIM matrix of 2nd order partial
!  derivatives.  Traverse the upper triangle of this matrix and,
!  for each element, compute a row of the least squares matrix.
!
do 123 idm = 1,mdim
    do 122 jdm = idm,mdim
        do 117 idim = 1,mdim
            nderiv(idim) = 0
117         continue
!
!  Off-diagonal elements appear twice by symmetry, so the corresponding
!  row is weighted by a factor of 2.
!
        rowwt = 2.d0*dcwght
        if (jdm/=idm) go to 118
!
!  Diagonal.
!
        rowwt = dcwght
        nderiv(jdm) = 2
        if (in(idm)/=0 .and. in(idm)/=inmx(idm)) go to 119
!
!  Node is at boundary.
!
!  Normal 2nd derivative constraint at boundary is not appropriate for
!  natural splines (2nd derivative 0 by definition).  Substitute
!  a 1st derivative constraint.
!
118         nderiv(idm) = 1
        nderiv(jdm) = 1
119         irow = irow + 1
!
!  Begining of basis index loop - traverse all indices corresponding
!  to basis functions which are non-zero at X.
!  The indices are in IB and are passed through common to BASCMP.
!
120         call bascmpd(x,nderiv,xmin,nodes,icol,basm)
!
!  BASCMP computes ICOL and BASM where BASM is the value at X of the
!  N-dimensional basis function corresponding to column ICOL.
!
        coef(icol) = rowwt*basm
!
!  Increment the basis indices.
!
        do 121 idim = 1,mdim
            ib(idim) = ib(idim) + 1
            if (ib(idim)<=ibmx(idim)) go to 120
            ib(idim) = ibmn(idim)
121         continue
!
!  End of basis index loop.
!
!  Send row of least squares matrix to reduction routine.
!
        call suprld(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef, &
                    reserr,lserr)
        if (lserr/=0) go to 133
122     continue
123 continue
!
!  Increment node indices.
!
124 do 125 idim = 1,mdim
    in(idim) = in(idim) + 1
    if (in(idim)<=inmx(idim)) go to 113
    in(idim) = 0
125 continue
!
!  End of node index loop.
!
!  Call for least squares solution in COEF array.
!
126 irow = 0
call suprld(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef,reserr, &
            lserr)
if (lserr/=0) go to 133
return
!
!  Error section
!
127 continue
ierror = 101
call cfaerr(ierror, &
    ' SPLCCD or SPLCWD - NDIM is less than 1 or is greater than 4' &
            ,60)
go to 134
128 continue
ierror = 102
call cfaerr(ierror, &
    ' SPLCCD or SPLCWD - NODES(IDIM) is less than 4 for some IDIM' &
            ,60)
go to 134
129 continue
ierror = 103
call cfaerr(ierror, &
  ' SPLCCD or SPLCWD - XMIN(IDIM) equals XMAX(IDIM) for some IDIM' &
            ,60)
go to 134
130 continue
ierror = 104
call cfaerr(ierror, &
    ' SPLCCD or SPLCWD - NCF (size of COEF) is too small         ' &
            ,60)
go to 134
131 continue
ierror = 105
call cfaerr(ierror, &
    ' SPLCCD or SPLCWD - Ndata Is less than 1                    ' &
            ,60)
go to 134
132 continue
ierror = 106
call cfaerr(ierror, &
    ' SPLCCD or SPLCWD - NWRK (size of WORK) is too small        ' &
            ,60)
go to 134
133 continue
ierror = 107
call cfaerr(ierror, &
' SPLCCD or SPLCWD - SUPRLS failure (this usually indicates insuff &
icient input data',80)
!
134 return
end
!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!
function splded(ndim,x,nderiv,coef,xmin,xmax,nodes,ierror)
double precision splded
double precision x
double precision coef
double precision xmin
double precision xmax
double precision dx
double precision dxin
double precision xrng
double precision sum
double precision basm
dimension x(ndim),nderiv(ndim),coef(*),xmin(ndim),xmax(ndim), &
          nodes(ndim)
common /splcomd/dx(4),dxin(4),mdim,ib(4),ibmn(4),ibmx(4)
save
!
! The restriction for NDIM to be <= 4 can be eliminated by increasing
! the above dimensions.
!
ierror = 0
mdim = ndim
if (mdim<1 .or. mdim>4) go to 105
iibmx = 1
do 101 idim = 1,mdim
    nod = nodes(idim)
    if (nod<4) go to 106
    xrng = xmax(idim) - xmin(idim)
    if (xrng==0.d0) go to 107
    if (nderiv(idim)<0 .or. nderiv(idim)>2) go to 108
!
!  DX(IDIM) is the node spacing along the IDIM coordinate.
!
    dx(idim) = xrng/dble(nod-1)
    dxin(idim) = 1.d0/dx(idim)
!
!  Compute indices of basis functions which are nonzero at X.
!
    it = dxin(idim)* (x(idim)-xmin(idim))
!
!  IBMN must be in the range 0 to NODES-2.
!
    ibmn(idim) = min0(max0(it-1,0),nod-2)
!
!  IBMX must be in the range 1 to NODES-1.
!
    ibmx(idim) = max0(min0(it+2,nod-1),1)
    iibmx = iibmx* (ibmx(idim)-ibmn(idim)+1)
    ib(idim) = ibmn(idim)
101 continue
!
sum = 0.d0
iib = 0
!
!  Begining of basis index loop - traverse all indices corresponding
!  to basis functions which are nonzero at X.
!
102 iib = iib + 1
!
!  The indices are in IB and are passed through common to BASCMP.
!
call bascmpd(x,nderiv,xmin,nodes,icof,basm)
!
!  BASCMP computes ICOF and BASM where BASM is the value at X of the
!  N-dimensional basis function corresponding to COEF(ICOF).
!
sum = sum + coef(icof)*basm
if (iib>=iibmx) go to 104
!
!  Increment the basis indices.
!
do 103 idim = 1,mdim
    ib(idim) = ib(idim) + 1
    if (ib(idim)<=ibmx(idim)) go to 102
    ib(idim) = ibmn(idim)
103 continue
!
!  End of basis index loop.
!
104 splded = sum
return
!
!  Errors.
!
105 continue
ierror = 101
call cfaerr(ierror, &
    ' SPLFED or SPLDED - NDIM is less than 1 or greater than 4   ' &
            ,60)
go to 109
106 continue
ierror = 102
call cfaerr(ierror, &
    ' SPLFED or SPLDED - NODES(IDIM) is less than  4for some IDIM' &
            ,60)
go to 109
107 continue
ierror = 103
call cfaerr(ierror, &
    ' SPLFED or SPLDED - XMIN(IDIM) = XMAX(IDIM) for some IDIM   ' &
            ,60)
go to 109
108 continue
ierror = 104
call cfaerr(ierror, &
' SPLDED - NDERIV(IDIM) IS less than 0 or greater than 2 for some &
IDIM  ',70)
!
109 stop
end
!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!
function splfed(ndim,x,coef,xmin,xmax,nodes,ierror)
double precision splfed
double precision x
double precision coef
double precision xmin
double precision xmax
double precision splded
dimension x(ndim),coef(*),xmin(ndim),xmax(ndim),nodes(ndim)
dimension nderiv(4)
save
!
data nderiv(1),nderiv(2),nderiv(3),nderiv(4)/0,0,0,0/
!
!  The restriction for NDIM to be <= 4 can be eliminated by
!  increasing the above dimension and those in SPLDED.
!
splfed = splded(ndim,x,nderiv,coef,xmin,xmax,nodes,ierror)
!
return
end
!
!                Copyright (C)  2000
!        University Corporation for Atmospheric Research
!                All Rights Reserved
!
! The use of this Software is governed by a License Agreement.
!
subroutine suprld(i,rowi,n,bi,a,nn,soln,err,ier)
double precision rowi
double precision bi
double precision a
double precision soln
double precision err
double precision errsum
double precision s
double precision temp
double precision temp1
double precision cn
double precision sn
dimension rowi(n),a(nn),soln(n)
save
!
ier = 0
if (i>1) go to 101
!
!  Routine entered with I<=0 means complete the reduction and store
!  the solution in SOLN.
!
if (i<=0) go to 125
!
!  Set up quantities on first call.
!
iold = 0
np1 = n + 1
!
!  Compute how many rows can be input now.
!
l = nn/np1
ilast = 0
il1 = 0
k = 0
k1 = 0
errsum = 0.d0
nreq = ((n+5)*n+2)/2
!
!  Error exit if insufficient scratch storage provided.
!
if (nn>=nreq) go to 101
ier = 32
call cfaerr(ier, &
' SUPRLD - insufficient scratch storage provided. at least ((N+5)* &
N+2)/2 locations needed',88)
return
!
!  Store the row in the scratch storage.
!
101 continue
!
!  Error exit if (I-IOLD)/=1.
!
if ((i-iold)==1) go to 102
ier = 35
call cfaerr(ier,' SUPRLD - values of I not in sequence',37)
return
!
102 continue
iold = i
do 103 j = 1,n
    ilj = ilast + j
    a(ilj) = rowi(j)
103 continue
ilnp = ilast + np1
a(ilnp) = bi
ilast = ilast + np1
isav = i
if (i<l) return
104 continue
if (k==0) go to 115
k1 = min0(k,n)
idiag = -np1
if (l-k==1) go to 110
!
!  Apply householder transformations to zero out new rows.
!
do 109 j = 1,k1
    idiag = idiag + (np1-j+2)
    i1 = il1 + j
    i2 = i1 + np1* (l-k-1)
    s = a(idiag)*a(idiag)
    do 105 ii = i1,i2,np1
        s = s + a(ii)*a(ii)
105     continue
    if (s==0.d0) go to 109
    temp = a(idiag)
    a(idiag) = sqrt(s)
    if (temp>0.d0) a(idiag) = -a(idiag)
    temp = temp - a(idiag)
    temp1 = 1.d0/ (temp*a(idiag))
    jp1 = j + 1
    do 108 j1 = jp1,np1
        jdel = j1 - j
        idj = idiag + jdel
        s = temp*a(idj)
        do 106 ii = i1,i2,np1
            iijd = ii + jdel
            s = s + a(ii)*a(iijd)
106         continue
        s = s*temp1
        a(idj) = a(idj) + s*temp
        do 107 ii = i1,i2,np1
            iijd = ii + jdel
            a(iijd) = a(iijd) + s*a(ii)
107         continue
108     continue
109 continue
go to 113
!
!  Apply rotations to zero out the single new row.
!
110 do 112 j = 1,k1
    idiag = idiag + (np1-j+2)
    i1 = il1 + j
    if (abs(a(i1))<=1.d-18) then
        s = sqrt(a(idiag)*a(idiag))
    else if (abs(a(idiag))<1.d-18) then
        s = sqrt(a(i1)*a(i1))
    else
        s = sqrt(a(idiag)*a(idiag)+a(i1)*a(i1))
    end if
    if (s==0.d0) go to 112
    temp = a(idiag)
    a(idiag) = s
    s = 1.d0/s
    cn = temp*s
    sn = a(i1)*s
    jp1 = j + 1
    do 111 j1 = jp1,np1
        jdel = j1 - j
        idj = idiag + jdel
        temp = a(idj)
        i1jd = i1 + jdel
        a(idj) = cn*temp + sn*a(i1jd)
        a(i1jd) = -sn*temp + cn*a(i1jd)
111     continue
112 continue
113 if (k<n) go to 115
lmkm1 = l - k
!
!  Accumulate residual sum of squares.
!
do 114 ii = 1,lmkm1
    ilnp = il1 + ii*np1
    errsum = errsum + a(ilnp)*a(ilnp)
114 continue
if (i<=0) go to 127
k = l
ilast = il1
!
!  Determine how many new rows may be input on next iteration.
!
l = k + (nn-ilast)/np1
return
115 k11 = k1 + 1
k1 = min0(l,n)
if (l-k==1) go to 122
k1m1 = k1 - 1
if (l>n) k1m1 = n
i1 = il1 + k11 - np1 - 1
!
!  Perform householder transformations to reduce rows to upper
!  triangular form.
!
do 120 j = k11,k1m1
    i1 = i1 + (np1+1)
    i2 = i1 + (l-j)*np1
    s = 0.d0
    do 116 ii = i1,i2,np1
        s = s + a(ii)*a(ii)
116     continue
    if (s==0.d0) go to 120
    temp = a(i1)
    a(i1) = sqrt(s)
    if (temp>0.d0) a(i1) = -a(i1)
    temp = temp - a(i1)
    temp1 = 1.d0/ (temp*a(i1))
    jp1 = j + 1
    i11 = i1 + np1
    do 119 j1 = jp1,np1
        jdel = j1 - j
        i1jd = i1 + jdel
        s = temp*a(i1jd)
        do 117 ii = i11,i2,np1
            iijd = ii + jdel
            s = s + a(ii)*a(iijd)
117         continue
        s = s*temp1
        i1jd = i1 + jdel
        a(i1jd) = a(i1jd) + s*temp
        do 118 ii = i11,i2,np1
            iijd = ii + jdel
            a(iijd) = a(iijd) + s*a(ii)
118         continue
119     continue
120 continue
if (l<=n) go to 122
np1mk = np1 - k
lmk = l - k
!
!  Accumulate residual sum of squares.
!
do 121 ii = np1mk,lmk
    ilnp = il1 + ii*np1
    errsum = errsum + a(ilnp)*a(ilnp)
121 continue
122 imov = 0
i1 = il1 + k11 - np1 - 1
!
!  Squeeze the unnecessary elements out of scratch storage to
!  allow space for more rows.
!
do 124 ii = k11,k1
    imov = imov + (ii-1)
    i1 = i1 + np1 + 1
    i2 = i1 + np1 - ii
    do 123 iii = i1,i2
        iiim = iii - imov
        a(iiim) = a(iii)
123     continue
124 continue
ilast = i2 - imov
il1 = ilast
if (i<=0) go to 127
k = l
!
!  Determine how many new rows may be input on next iteration.
!
l = k + (nn-ilast)/np1
return
!
!  Complete reduction and store solution in SOLN.
!
125 l = isav
!
!  Error exit if L less than N.
!
if (l>=n) go to 126
ier = 33
call cfaerr(ier,' SUPRLD - array has too few rows.',33)
return
126 continue
!
!  K/=ISAV means further reduction needed.
!
if (k/=isav) go to 104
127 ilast = (np1* (np1+1))/2 - 1
if (a(ilast-1)==0.d0) go to 130
!
! Solve triangular system into ROWI.
!
soln(n) = a(ilast)/a(ilast-1)
do 129 ii = 2,n
    iim1 = ii - 1
    ilast = ilast - ii
    s = a(ilast)
    do 128 k = 1,iim1
        ilk = ilast - k
        npk = np1 - k
        s = s - a(ilk)*soln(npk)
128     continue
    ilii = ilast - ii
    if (a(ilii)==0.d0) go to 130
    npii = np1 - ii
    soln(npii) = s/a(ilii)
129 continue
!
!  Store residual norm.
!
err = sqrt(errsum)
return
!
!  Error return if system is singular.
!
130 continue
ier = 34
call cfaerr(ier,' SUPRLD - system is singular.',29)
return
!
end

    end module splpak_module
