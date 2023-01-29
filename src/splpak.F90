!*****************************************************************************************
!>
!  This package contains routines for fitting
!  (least squares) a multidimensional cubic spline
!  to arbitrarily located data.  It also contains
!  routines for evaluating this spline (or its
!  partial derivatives) at any point.
!
!  Coefficient calculation is performed in
!  subroutines [[splcc]] or [[splcw]] and evaluation is
!  performed by functions [[splfe]] or [[splde]].
!
!### History
!  * Developed in 1972-73 by Dave Fulker of NCAR's Scientific Computing Division.
!  * Cleaned up and added to the Ngmath library in 1998.
!  * Latest revision to the original Fortran code: August, 1998
!  * Jacob Williams : Jan 2023 : modernized this code.
!
!### License
!    Copyright (C) 2000
!    University Corporation for Atmospheric Research
!    All Rights Reserved
!    The use of this Software is governed by a License Agreement.

    module splpak_module

    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter :: wp = real32   !! Real working precision [4 bytes]
#elif REAL64
    integer,parameter :: wp = real64   !! Real working precision [8 bytes]
#elif REAL128
    integer,parameter :: wp = real128  !! Real working precision [16 bytes]
#else
    integer,parameter :: wp = real64   !! Real working precision if not specified [8 bytes]
#endif

    integer,parameter,public :: splpak_wp = wp   !! Working precision

    type,public :: splpak_type

        !!### Usage
        !!
        !!  The class contains four user entries:
        !!  [[splcc]], [[splcw]], [[splfe]], and [[splde]].
        !!
        !!  The user first calls [[splcc]] by
        !!```fortran
        !!    call me%initialize(ndim,xdata,l1xdat,ydata,ndata,
        !!                       xmin,xmax,nodes,xtrap,coef,ncf,
        !!                       work,nwrk,ierror)
        !!```
        !!  or [[splcw]] by
        !!```fortran
        !!    call me%initialize(ndim,xdata,l1xdata,ydata,wdata,
        !!                       ndata,xmin,xmax,nodes,xtrap,
        !!                       coef,ncf,work,nwrk,ierror)
        !!```
        !!  The parameter `NDATA` in the call to [[splcw]]
        !!  enables the user to weight some of the data
        !!  points more heavily than others.  Both
        !!  routines return a set of coefficients in the
        !!  array `COEF`.  These coefficients are
        !!  subsequently used in the computation of
        !!  function values and partial derivatives.
        !!  To compute values on the spline approximation
        !!  the user then calls [[splfe]] or [[splde]] any
        !!  number of times in any order provided that
        !!  the values of the inputs, `NDIM`, `COEF`, `XMIN`,
        !!  `XMAX`, and `NODES`, are preserved between calls.
        !!
        !!  [[splfe]] and [[splde]] are called in the following way:
        !!```fortran
        !!    f = me%evaluate(ndim,x,coef,xmin,xmax,nodes,ierror)
        !!```
        !!  or
        !!```fortran
        !!    f = me%evaluate(ndim,x,nderiv,coef,xmin,xmax,nodes,ierror)
        !!```
        !!  The routine [[splfe]] returns an interpolated
        !!  value at the point defined by the array `X`.
        !!  [[splde]] affords the user the additional
        !!  capability of calculating an interpolated
        !!  value for one of several partial derivatives
        !!  specified by the array `NDERIV`.

        private

        ! formerly in splcomd common block:
        integer :: mdim = 0
        real(wp),dimension(:),allocatable :: dx ! originally these were all size 4
        real(wp),dimension(:),allocatable :: dxin
        integer,dimension(:),allocatable  :: ib
        integer,dimension(:),allocatable  :: ibmn
        integer,dimension(:),allocatable  :: ibmx

        ! formerly saved variables in suprls
        integer :: ilast = 0
        integer :: isav  = 0
        integer :: iold  = 0
        integer :: np1   = 0
        integer :: l     = 0
        integer :: il1   = 0
        integer :: k     = 0
        integer :: k1    = 0
        real(wp) :: errsum = 0.0_wp

        contains

        private

        generic,public   :: initialize => splcc, splcw !! compute the spline coefficients
        generic,public   :: evaluate   => splfe, splde !! evaluate the spline
        procedure,public :: destroy    => destroy_splpak !! destory the internal class variables
        procedure,private :: splcc
        procedure,private :: splcw
        procedure,private :: splfe
        procedure,private :: splde
        procedure,private :: bascmp
        procedure,private :: suprls

    end type splpak_type

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Destroy the internal class variables.

    subroutine destroy_splpak(me,ndim)

        class(splpak_type),intent(inout) :: me
        integer,intent(in),optional :: ndim

        if (allocated(me%dx))   deallocate(me%dx)
        if (allocated(me%dxin)) deallocate(me%dxin)
        if (allocated(me%ib))   deallocate(me%ib)
        if (allocated(me%ibmn)) deallocate(me%ibmn)
        if (allocated(me%ibmx)) deallocate(me%ibmx)
        if (present(ndim)) then
            allocate(me%dx  (ndim)); me%dx   = 0.0_wp
            allocate(me%dxin(ndim)); me%dxin = 0.0_wp
            allocate(me%ib  (ndim)); me%ib   = 0
            allocate(me%ibmn(ndim)); me%ibmn = 0
            allocate(me%ibmx(ndim)); me%ibmx = 0
        end if

        me%mdim   = 0
        me%ilast  = 0
        me%isav   = 0
        me%iold   = 0
        me%np1    = 0
        me%l      = 0
        me%il1    = 0
        me%k      = 0
        me%k1     = 0
        me%errsum = 0.0_wp

    end subroutine destroy_splpak
!*****************************************************************************************

!*****************************************************************************************
!>
!  This routine does basis function computations for natural
!  splines.  This routine is called by routines [[SPLCW]] and [[SPLDE]]
!  to compute `ICOL` and `BASM`, which are defined as follows:
!
!   * The `MDIM` indices in `IB` (defined in [[splpak_type]]) determine
!     a specific node in the node grid (see routine [[SPLCC]] for a
!     description of the node grid).  Every node is associated
!     with an `MDIM`-dimensional basis function and a corresponding
!     column in the least squares matrix (or element of the
!     coefficient vector).  The column index (which may be thought
!     of as a linear address for the `MDIM`-dimensional node grid)
!     corresponding to the specified node is computed as `ICOL`.  The
!     associated basis function evaluated at `X` (an `MDIM`-vector) is
!     computed as `BASM` (a scalar).
!
!  In case `NDERIV` is not all zero, `BASM` will not be the value of
!  the basis function but rather a partial derivative of that
!  function as follows:
!
!   * The order of the partial derivative in the direction of the
!     `IDIM` coordinate is `NDERIV(IDIM)` (for `IDIM <= MDIM`).  This
!     routine will compute incorrect values if `NDERIV(IDIM)` is not
!     in the range 0 to 2.
!
!  The technique of this routine is to transform the independent
!  variable in each dimension such that the nodes fall on
!  suitably chosen integers.  On this transformed space, the
!  1-dimensional basis functions and their derivatives have a
!  particularly simple form.  The desired `MDIM`-dimensional basis
!  function (or any of its partial derivatives) is computed as
!  a product of such 1-dimensional functions (tensor product
!  method of defining multi-dimensional splines).  The values
!  which determine the location of the nodes, and hence the
!  above transform, are passed through common and the argument
!  list.

subroutine bascmp(me,x,nderiv,xmin,nodes,icol,basm)

    class(splpak_type),intent(inout) :: me
    real(wp),intent(in) :: x(:)
    integer,intent(in) :: nderiv(:)
    real(wp),intent(in) :: xmin(:)
    integer,intent(in) :: nodes(:)
    integer,intent(out) :: icol
    real(wp),intent(out) :: basm

    real(wp) :: xb,bas1,z,fact,z1
    integer :: idim,mdmid,ntyp,ngo

    ! ICOL will be a linear address corresponding to the indices in IB.
    icol = 0

    ! BASM will be M-dimensional basis function evaluated at X.
    basm = 1.0_wp
    do idim = 1,me%mdim

        ! Compute ICOL by Horner's method.
        mdmid = me%mdim + 1 - idim
        icol = nodes(mdmid)*icol + me%ib(mdmid)

        ! NGO depends upon function type and NDERIV.
        ntyp = 1

        ! Function type 1 (left linear) for IB = 0 or 1.
        if (me%ib(idim)>1) then
            ntyp = 2
            ! Function type 2 (chapeau function) for 2 LT IB LT NODES-2.
            if (me%ib(idim)>=nodes(idim)-2) then
                ntyp = 3
            end if
        end if

        ! Function type 3 (right linear) for IB = NODES-2 or NODES-1.
        ngo = 3*ntyp + nderiv(idim) - 2

        !  XB is X value of node IB (center of basis function).
        xb = xmin(idim) + real(me%ib(idim),wp)*me%dx(idim)

        !  BAS1 will be the 1-dimensional basis function evaluated at X.
        bas1 = 0.0_wp

        select case (ngo)

        case(4)

            !  Function type 2 (chapeau function).
            !
            !  Transform so that XB is at the origin and the other nodes are at
            !  the integers.
            z = abs(me%dxin(idim)* (x(idim)-xb)) - 2.0_wp

            !  This chapeau function is then that unique cubic spline which is
            !  identically zero for ABS(Z) GE 2 and is 1 at the origin.  This
            !  function is the general interior node basis function.
            if (z<0.0_wp) then
                bas1 = -0.25_wp*z**3
                z = z + 1.0_wp
                if (z<0.0_wp) then
                    bas1 = bas1 + z**3
                end if
            end if

        case(5)

            !  1st derivative.
            z = x(idim) - xb
            fact = me%dxin(idim)
            if (z<0.0_wp) fact = -fact
            z = fact*z - 2.0_wp
            if (z<0.0_wp) then
                bas1 = -0.75_wp*z**2
                z = z + 1.0_wp
                if (z<0.0_wp) then
                    bas1 = bas1 + 3.0_wp*z**2
                end if
                bas1 = fact*bas1
            end if

        case(6) ! 108

            !  2nd derivative.
            fact = me%dxin(idim)
            z = fact*abs(x(idim)-xb) - 2.0_wp
            if (z<0.0_wp) then
                bas1 = -1.5_wp*z
                z = z + 1.0_wp
                if (z<0.0_wp) then
                    bas1 = bas1 + 6.0_wp*z
                end if
                bas1 = (fact**2)*bas1
            end if

        case(2,8)

            !  1st derivative.
            if (ngo==2) then
                fact = -me%dxin(idim)
            else if (ngo==8) then
                fact = me%dxin(idim)
            end if
            z = fact* (x(idim)-xb) + 2.0_wp
            if (z<0.0_wp) then
                if (z<2.0_wp) then
                    bas1 = 1.5_wp*z**2
                    z = z - 1.0_wp
                    if (z>0.0_wp) then
                        bas1 = bas1 - 3.0_wp*z**2
                    end if
                    bas1 = fact*bas1
                else
                    bas1 = 3.0_wp*fact
                end if
            end if

        case(3,9)

            !  2nd derivative.
            if (ngo==3) then
                fact = -me%dxin(idim)
            else if (ngo==9) then
                fact = me%dxin(idim)
            end if
            z = fact* (x(idim)-xb) + 2.0_wp
            z1 = z - 1.0_wp
            if (abs(z1)<1.0_wp) then
                bas1 = 3.0_wp*z
                if (z1>0.0_wp) then
                    bas1 = bas1 - 6.0_wp*z1
                end if
                bas1 = (fact**2)*bas1
            end if

        case default ! case(1,7) ! or ngo some other value (does that ever happen?)
                     !             (due to the computed goto in the original code)

            if (ngo/=7) then
                !  Function type 1 (left linear) is mirror image of function type 3.
                !
                !  Transform so that XB is at 2 and the other nodes are at the integers
                !  (with ordering reversed to form a mirror image).
                z = me%dxin(idim)* (xb-x(idim)) + 2.0_wp
            else
                !  Function type 3 (right linear).
                !
                !  Transform so that XB is at 2 and the other nodes are at the integers.
                z = me%dxin(idim)* (x(idim)-xb) + 2.0_wp
            end if

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

            if (z>0.0_wp) then
                if (z<2.0_wp) then
                    bas1 = 0.5_wp*z**3
                    z = z - 1.0_wp
                    if (z>0.0_wp) then
                        bas1 = bas1 - z**3
                    end if
                else
                    bas1 = 3.0_wp*z - 3.0_wp
                end if
            end if

        end select

        basm = basm*bas1

    end do

    icol = icol + 1

end subroutine bascmp
!*****************************************************************************************

!*****************************************************************************************
!>
!  To print an error number and an error message
!  or just an error message.
!
!  The message is writen to `output_unit`.

subroutine cfaerr (ierr,mess)

    integer,intent(in) :: ierr !! The error number (printed only if non-zero).
    character(len=*), intent(in) :: mess !! Message to be printed.

    if (ierr /= 0) write (output_unit,'(A,I5)') ' IERR=', ierr
    write (output_unit,'(A)') trim(mess)

end subroutine cfaerr
!*****************************************************************************************

!*****************************************************************************************
!>
!  N-dimensional cubic spline coefficient
!  calculation by least squares.
!
!  The usage and arguments of this routine are
!  identical to those for [[SPLCW]] except for the
!  omission of the array of weights, `WDATA`.  See
!  entry [[SPLCW]] description for a
!  complete description.

subroutine splcc(me,ndim,xdata,l1xdat,ydata,ndata,xmin,xmax,nodes, &
                 xtrap,coef,ncf,work,nwrk,ierror)

    class(splpak_type),intent(inout) :: me
    integer,intent(in) :: ndim
    integer,intent(in) :: l1xdat
    integer,intent(in) :: ncf
    integer,intent(in) :: nwrk
    integer,intent(in) :: ndata
    real(wp),intent(in) :: xdata(l1xdat,ndata)
    real(wp),intent(in) :: ydata(ndata)
    real(wp),intent(in) :: xmin(ndim)
    real(wp),intent(in) :: xmax(ndim)
    real(wp),intent(in) :: xtrap
    integer,intent(in) :: nodes(ndim)
    real(wp) :: work(nwrk)
    real(wp),intent(out) :: coef(ncf)
    integer,intent(out) :: ierror

    real(wp),dimension(1),parameter :: wdata = -1.0_wp !! indicates to [[splcw]]
                                                       !! that weights are not used

    call me%splcw(ndim,xdata,l1xdat,ydata,wdata,ndata,xmin,xmax,&
                  nodes,xtrap,coef,ncf,work,nwrk,ierror)

end subroutine splcc
!*****************************************************************************************

!*****************************************************************************************
!>
!  N-dimensional cubic spline coefficient
!  calculation by weighted least squares on
!  arbitrarily located data.
!
!  The spline (or its derivatives) may then be
!  evaluated by using function [[SPLFE]] (or [[SPLDE]]).
!
!  A grid of evenly spaced nodes in NDIM space is
!  defined by the arguments XMIN, XMAX and NODES.
!  A linear basis for the class of natural splines
!  on these nodes is formed, and a set of
!  corresponding coefficients is computed in the
!  array COEF.  These coefficients are chosen to
!  minimize the weighted sum of squared errors
!  between the spline and the arbitrarily located
!  data values described by the arguments XDATA,
!  YDATA and NDATA.  The smoothness of the spline
!  in data sparse areas is controlled by the
!  argument XTRAP.
!
!### Note
!  In order to understand the arguments of this
!  routine, one should realize that the node grid
!  need not bear any particular relation to the
!  data points.  In the theory of exact-fit
!  interpolatory splines, the nodes would in fact
!  be data locations, but in this case they serve
!  only to define the class of splines from which
!  the approximating function is chosen.  This
!  node grid is a rectangular arrangement of
!  points in NDIM space, with the restriction that
!  along any coordinate direction the nodes are
!  equally spaced.  The class of natural splines
!  on this grid of nodes (NDIM-cubic splines whose
!  2nd derivatives normal to the boundaries are 0)
!  has as many degrees of freedom as the grid has
!  nodes.  Thus the smoothness or flexibility of
!  the splines is determined by the choice of the
!  node grid.
!
!### Algorithm
!  An overdetermined system of linear equations
!  is formed -- one equation for each data point
!  plus equations for derivative constraints.
!  This system is solved using subroutine [[suprls]].
!
!### Accuracy
!  If there is exactly one data point in the
!  near vicinity of each node and no extra data,
!  the resulting spline will agree with the
!  data values to machine accuracy.  However, if
!  the problem is overdetermined or the sparse
!  data option is utilized, the accuracy is hard
!  to predict.  Basically, smooth functions
!  require fewer nodes than rough ones for the
!  same accuracy.
!
!### Timing
!  The execution time is roughly proportional
!  to `NDATA*NCOF**2` where `NCOF = NODES(1)*...*NODES(NDIM)`.

subroutine splcw(me,ndim,xdata,l1xdat,ydata,wdata,ndata,xmin,xmax, &
                 nodes,xtrap,coef,ncf,work,nwrk,ierror)

    class(splpak_type),intent(inout) :: me
    integer,intent(in) :: ndim !! The dimensionality of the problem.  The
                               !! spline is a function of `NDIM` variables or
                               !! coordinates and thus a point in the
                               !! independent variable space is an `NDIM` vector.
                               !! `NDIM` must be `>= 1`.
    integer,intent(in) :: l1xdat !! The length of the 1st dimension of `XDATA` in
                                 !! the calling program.  `L1XDAT` must be `>= NDIM`.
                                 !!
                                 !!#### Note:
                                 !! For 1-dimensional problems `L1XDAT` is usually 1.
    integer,intent(in) :: ncf !! The length of the array `COEF` in the calling
                              !! program.  If `NCF` is `< NODES(1)*...*NODES(NDIM)`,
                              !! a fatal error is diagnosed.
    integer,intent(in) :: nwrk !! The length of the array `WORK` in the calling
                               !! program.  If
                               !! `NCOL = NODES(1)*...*NODES(NDIM)` is the total
                               !! number of nodes, then a fatal error is
                               !! diagnosed if `NWRK` is less than
                               !! `NCOL*(NCOL+1)`.
    integer,intent(in) :: ndata !! The number of data points mentioned in the
                                !! above arguments.
    real(wp),intent(in) :: xdata(l1xdat,ndata) !! A collection of locations for the data
                                               !! values, i.e., points from the independent
                                               !! variable space.  This collection is a
                                               !! 2-dimensional array whose 1st dimension
                                               !! indexes the `NDIM` coordinates of a given point
                                               !! and whose 2nd dimension labels the data
                                               !! point.  For example, the data point with
                                               !! label `IDATA` is located at the point
                                               !! `(XDATA(1,IDATA),...,XDATA(NDIM,IDATA))` where
                                               !! the elements of this vector are the values of
                                               !! the `NDIM` coordinates.  The location, number
                                               !! and ordering of the data points is arbitrary.
                                               !! The dimension of `XDATA` is assumed to be
                                               !! `XDATA(L1XDAT,NDATA)`.
    real(wp),intent(in) :: ydata(ndata) !! A collection of data values corresponding to
                                        !! the points in `XDATA`.  `YDATA(IDATA)` is the
                                        !! data value associated with the point
                                        !! `(XDATA(1,IDATA),...,XDATA(NDIM,IDATA))` in the
                                        !! independent variable space.  The spline whose
                                        !! coefficients are computed by this routine
                                        !! approximates these data values in the least
                                        !! squares sense.  The dimension is assumed to be
                                        !! `YDATA(NDATA)`.
    real(wp),intent(in) :: wdata(:) !! A collection of weights.  `WDATA(IDATA)` is a
                                    !! weight associated with the data point
                                    !! labelled `IDATA`.  It should be non-negative,
                                    !! but may be of any magnitude.  The weights
                                    !! have the effect of forcing greater or lesser
                                    !! accuracy at a given point as follows: this
                                    !! routine chooses coefficients to minimize the
                                    !! sum over all data points of the quantity
                                    !!```fortran
                                    !!   (WDATA(IDATA)*(YDATA(IDATA) ! spline value at XDATA(IDATA)))**2.
                                    !!```
                                    !! Thus, if the reliability
                                    !! of a data point is known to be low, the
                                    !! corresponding weight may be made small
                                    !! (relative to the other weights) so that the
                                    !! sum over all data points is affected less by
                                    !! discrepencies at the unreliable point.  Data
                                    !! points with zero weight are completely
                                    !! ignored.
                                    !!
                                    !!#### Note:
                                    !! If `WDATA(1)` is `< 0`, the other
                                    !! elements of `WDATA` are not
                                    !! referenced, and all weights are
                                    !! assumed to be unity.
                                    !!
                                    !! The dimension is assumed to be `WDATA(NDATA)`
                                    !! unless `WDATA(1) < 0.`, in which case the
                                    !! dimension is assumed to be 1.
    real(wp),intent(in) :: xmin(ndim) !! A vector describing the lower extreme corner
                                      !! of the node grid.  A set of evenly spaced
                                      !! nodes is formed along each coordinate axis
                                      !! and `XMIN(IDIM)` is the location of the first
                                      !! node along the `IDIM` axis.  The dimension is
                                      !! assumed to be `XMIN(NDIM)`.
    real(wp),intent(in) :: xmax(ndim) !! A vector describing the upper extreme corner
                                      !! of the node grid.  A set of evenly spaced
                                      !! nodes is formed along each coordinate axis
                                      !! and `XMAX(IDIM)` is the location of the last
                                      !! node along the `IDIM` axis.  The dimension is
                                      !! assumed to be `XMAX(NDIM)`.
    real(wp),intent(in) :: xtrap !! A parameter to control extrapolation to data
                                 !! sparse areas.  The region described by `XMIN`
                                 !! and `XMAX` is divided into rectangles, the
                                 !! number of which is determined by `NODES`, and
                                 !! any rectangle containing a disproportionately
                                 !! small number of data points is considered to
                                 !! be data sparse (rectangle is used here to
                                 !! mean `NDIM`-dimensional rectangle).  If `XTRAP`
                                 !! is nonzero the least squares problem is
                                 !! augmented with derivative constraints in the
                                 !! data sparse areas to prevent the matrix from
                                 !! becoming poorly conditioned.  `XTRAP` serves as
                                 !! a weight for these constraints, and thus may
                                 !! be used to control smoothness in data sparse
                                 !! areas.  Experience indicates that unity is a
                                 !! good first guess for this parameter.
                                 !!
                                 !!#### Note:
                                 !! If `XTRAP` is zero, substantial
                                 !! portions of the routine will be
                                 !! skipped, but a singular matrix
                                 !! can result if large portions of
                                 !! the region are without data.
    integer,intent(in) :: nodes(ndim) !! A vector of integers describing the number of
                                      !! nodes along each axis.  `NODES(IDIM)` is the
                                      !! number of nodes (counting endpoints) along
                                      !! the `IDIM` axis and determines the flexibility
                                      !! of the spline in that coordinate direction.
                                      !! `NODES(IDIM)` must be `>= 4`, but may be as
                                      !! large as the arrays `COEF` and `WORK` allow.
                                      !! The dimension is assumed to be `NODES(NDIM)`.
                                      !!
                                      !!#### Note:
                                      !! The node grid is completely defined by
                                      !! the arguments `XMIN`, `XMAX` and `NODES`.
                                      !! The spacing of this grid in the `IDIM`
                                      !! coordinate direction is:
                                      !!```fortran
                                      !!   DX(IDIM) = (XMAX(IDIM)-XMIN(IDIM)) / (NODES(IDIM)-1).
                                      !!```
                                      !! A node in this grid may be indexed by
                                      !! an `NDIM` vector of integers
                                      !! `(IN(1),...,IN(NDIM))` where
                                      !! `1 <= IN(IDIM) <= NODES(IDIM)`.
                                      !! The location of such a node may be
                                      !! represented by an `NDIM` vector
                                      !! `(X(1),...,X(NDIM))` where
                                      !! `X(IDIM) = XMIN(IDIM) + (IN(IDIM)-1) * DX(IDIM)`.
    real(wp) :: work(nwrk) !! A workspace array for solving the least
                           !! squares matrix generated by this routine.
                           !! Its required size is a function of the total
                           !! number of nodes in the node grid.  This
                           !! total, `NCOL = NODES(1)*...*NODES(NDIM)`, is
                           !! also the number of columns in the least
                           !! squares matrix.  The length of the array `WORK`
                           !! must equal or exceed `NCOL*(NCOL+1)`.
    real(wp),intent(out) :: coef(ncf) !! The array of coefficients computed by this
                                      !! routine.  Each coefficient corresponds to a
                                      !! particular basis function which in turn
                                      !! corresponds to a node in the node grid.  This
                                      !! correspondence between the node grid and the
                                      !! array `COEF` is as if `COEF` were an
                                      !! `NDIM`-dimensional Fortran array with
                                      !! dimensions `NODES(1),...,NODES(NDIM)`, i.e., to
                                      !! store the array linearly, the leftmost
                                      !! indices are incremented most frequently.
                                      !! Hence the length of the `COEF` array must equal
                                      !! or exceed the total number of nodes, which is
                                      !! `NODES(1)*...*NODES(NDIM)`.  The computed array
                                      !! `COEF` may be used with function [[SPLFE]]
                                      !! (or [[SPLDE]]) to evaluate the spline (or its
                                      !! derivatives) at an arbitrary point in `NDIM`
                                      !! space.  The dimension is assumed to be `COEF(NCF)`.
    integer,intent(out) :: ierror !! An error flag with the following meanings:
                                  !!
                                  !! * `  0`  No error.
                                  !! * `101`  `NDIM` is < 1.
                                  !! * `102`  `NODES(IDIM)` is < 4 for some `IDIM`.
                                  !! * `103`  `XMIN(IDIM) = XMAX(IDIM)` for some `IDIM`.
                                  !! * `104`  `NCF` (size of `COEF`) is `< NODES(1)*...*NODES(NDIM)`.
                                  !! * `105`  `NDATA` is `< 1`.
                                  !! * `106`  `NWRK` (size of `WORK`) is too small.
                                  !! * `107`  [[suprls]] failure (usually insufficient
                                  !!   data) -- ordinarily occurs only if
                                  !!   `XTRAP` is zero or `WDATA` contains all
                                  !!   zeros.

    real(wp),dimension(:),allocatable :: x
    integer,dimension(:),allocatable :: nderiv,in,inmx
    real(wp) :: xrng,swght,rowwt,rhs,basm,reserr,totlwt,&
                bump,wtprrc,expect,dcwght
    integer :: ncol,idim,nod,nwrk1,mdata,nwlft,irow,idata,&
               icol,it,lserr,iin,nrect,idimc,idm,jdm,inidim
    logical :: boundary

    real(wp),parameter :: spcrit = 0.75_wp
        !! SPCRIT is used to determine data sparseness as follows -
        !! the weights assigned to all data points are totaled into the
        !! variable TOTLWT. (If no weights are entered, it is set to
        !! NDATA.)  Each node of the node network is assigned a
        !! rectangle (in which it is contained) and the weights of all
        !! data points which fall in that rectangle are totaled.  If that
        !! total is less than SPCRIT*EXPECT (EXPECT is defined below),
        !! then the node is ascertained to be in a data sparse location.
        !! EXPECT is that fraction of TOTLWT that would be expected by
        !! comparing the area of the rectangle with the total area under
        !! consideration.

    ! size the arrays:
    call me%destroy(ndim)
    allocate(x(ndim))
    allocate(nderiv(ndim))
    allocate(in(ndim))
    allocate(inmx(ndim))

    ierror = 0
    me%mdim = ndim
    if (me%mdim<1) then
        ierror = 101
        call cfaerr(ierror, &
            ' splcc or splcw - NDIM is less than 1')
        return
    end if

    ncol = 1
    do idim = 1,me%mdim
        nod = nodes(idim)
        if (nod<4) then
            ierror = 102
            call cfaerr(ierror, &
                ' splcc or splcw - NODES(IDIM) is less than 4 for some IDIM')
            return
        end if

        !  Number of columns in least squares matrix = number of coefficients =
        !  product of nodes over all dimensions.
        ncol = ncol*nod
        xrng = xmax(idim) - xmin(idim)
        if (xrng==0.0_wp) then
            ierror = 103
            call cfaerr(ierror, &
                ' splcc or splcw - XMIN(IDIM) equals XMAX(IDIM) for some IDIM')
            return
        end if

        !  DX(IDIM) is the node spacing along the IDIM coordinate.
        me%dx(idim) = xrng/real(nod-1,wp)
        me%dxin(idim) = 1.0_wp/me%dx(idim)
        nderiv(idim) = 0
    end do
    if (ncol>ncf) then
        ierror = 104
        call cfaerr(ierror, &
            ' splcc or splcw - NCF (size of COEF) is too small')
        return
    end if
    nwrk1 = 1
    mdata = ndata
    if (mdata<1) then
        ierror = 105
        call cfaerr(ierror, &
            ' splcc or splcw - Ndata Is less than 1')
        return
    end if

    !  SWGHT is a local variable = XTRAP, and can be considered a smoothing
    !  weight for data sparse areas.  If SWGHT == 0, no smoothing
    !  computations are performed.
    swght = xtrap

    !  Set aside workspace for counting data points.
    if (swght/=0.0_wp) nwrk1 = ncol + 1

    !  NWLFT is the length of the remaining workspace.
    nwlft = nwrk - nwrk1 + 1
    if (nwlft<1) then
        ierror = 106
        call cfaerr(ierror, &
            ' splcc or splcw - NWRK (size of WORK) is too small')
        return
    end if
    irow = 0

    !  ROWWT is used to weight rows of the least squares matrix.
    rowwt = 1.0_wp

    !  Loop through all data points, computing a row for each.
    do idata = 1,mdata

        !  WDATA(1)<0 means weights have not been entered.  In that case,
        !  ROWWT is left equal to  1. for all points.  Otherwise ROWWT is
        !  equal to WDATA(IDATA).
        !
        !  Every element of the row, as well as the corresponding right hand
        !  side, is multiplied by ROWWT.
        if (wdata(1)>=0.0_wp) then
            rowwt = wdata(idata)
            !  Data points with 0 weight are ignored.
            if (rowwt==0.0_wp) cycle
        end if
        irow = irow + 1

        !  One row of the least squares matrix corresponds to each data
        !  point.  The right hand for that row will correspond to the
        !  function value YDATA at that point.
        rhs = rowwt*ydata(idata)
        do idim = 1,me%mdim
            x(idim) = xdata(idim,idata)
        end do

        !  The COEF array serves as a row of least squares matrix.
        !  Its value is zero except for columns corresponding to functions
        !  which are nonzero at X.
        do icol = 1,ncol
            coef(icol) = 0.0_wp
        end do

        !  Compute the indices of basis functions which are nonzero at X.
        !  IBMN is in the range 0 to nodes-2 and IBMX is in range 1
        !  to NODES-1.
        do idim = 1,me%mdim
            nod = nodes(idim)
            it = me%dxin(idim)* (x(idim)-xmin(idim))
            me%ibmn(idim) = min(max(it-1,0),nod-2)
            me%ib(idim) = me%ibmn(idim)
            me%ibmx(idim) = max(min(it+2,nod-1),1)
        end do

        basis_index : do
            !  Begining of basis index loop - traverse all indices corresponding
            !  to basis functions which are nonzero at X.  The indices are in
            !  IB and are passed through common to BASCMP.
            call me%bascmp(x,nderiv,xmin,nodes,icol,basm)

            !  BASCMP computes ICOL and BASM where BASM is the value at X of
            !  the N-dimensional basis function corresponding to column ICOL.
            coef(icol) = rowwt*basm

            !  Increment the basis indices.
            do idim = 1,me%mdim
                me%ib(idim) = me%ib(idim) + 1
                if (me%ib(idim)<=me%ibmx(idim)) cycle basis_index
                me%ib(idim) = me%ibmn(idim)
            end do
            exit basis_index !  End of basis index loop.
        end do basis_index

        !  Send a row of the least squares matrix to the reduction routine.
        call me%suprls(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef,reserr,lserr)
        if (lserr/=0) then
            ierror = 107
            call cfaerr(ierror, ' splcc or splcw - suprls failure '//&
                                '(this usually indicates insufficient input data)')
        end if
    end do

    !  Row computations for all data points are now complete.
    !
    !  If SWGHT==0, the least squares matrix is complete and no
    !  smoothing rows are computed.

    if (swght/=0.0_wp) then

        !  Initialize smoothing computations for data sparse areas.
        !  Derivative constraints will always have zero right hand side.
        rhs = 0.0_wp
        nrect = 1

        !  Initialize the node indices and compute number of rectangles
        !  formed by the node network.
        do idim = 1,me%mdim
            in(idim) = 0
            inmx(idim) = nodes(idim) - 1
            nrect = nrect*inmx(idim)
        end do

        !  Every node is assigned an element of the workspace (set aside
        !  previously) in which data points are counted.
        do iin = 1,ncol
            work(iin) = 0.0_wp
        end do

        !  Assign each data point to a node, total the assignments for
        !  each node, and save in the workspace.
        totlwt = 0.0_wp
        do idata = 1,mdata

            ! BUMP is the weight associated with the data point.
            bump = 1.0_wp
            if (wdata(1)>=0.0_wp) bump = wdata(idata)
            if (bump==0.0_wp) cycle

            ! Find the nearest node.
            iin = 0
            do idimc = 1,me%mdim
                idim = me%mdim + 1 - idimc
                inidim = int(me%dxin(idim)* (xdata(idim,idata)-xmin(idim))+0.5_wp)
                ! Points not in range (+ or - 1/2 node spacing) are not counted.
                if (inidim<0 .or. inidim>inmx(idim)) cycle
                ! Compute linear address of node in workspace by Horner's method.
                iin = (inmx(idim)+1)*iin + inidim
            end do

            ! Bump counter for that node.
            work(iin+1) = work(iin+1) + bump
            totlwt = totlwt + bump
        end do

        ! Compute the expected weight per rectangle.
        wtprrc = totlwt/real(nrect,wp)

        !  IN contains indices of the node (previously initialized).
        !  IIN will be the linear address of the node in the workspace.
        iin = 0

        !  Loop through all nodes, computing derivative constraint rows
        !  for those in data sparse locations.
        !
        !  Begining of node index loop - traverse all node indices.
        !  The indices are in IN.
        node_index : do
            iin = iin + 1
            expect = wtprrc

            !  Rectangles at edge of network are smaller and hence less weight
            !  should be expected.
            do idim = 1,me%mdim
                if (in(idim)==0 .or. in(idim)==inmx(idim)) expect = 0.5_wp*expect
            end do

            !  The expected weight minus the actual weight serves to define
            !  data sparseness and is also used to weight the derivative
            !  constraint rows.
            !
            !  There is no constraint if not data sparse.
            if (work(iin)<spcrit*expect) then

                dcwght = expect - work(iin)
                do idim = 1,me%mdim
                    inidim = in(idim)

                    !  Compute the location of the node.
                    x(idim) = xmin(idim) + real(inidim,wp)*me%dx(idim)

                    !  Compute the indices of the basis functions which are non-zero
                    !  at the node.
                    me%ibmn(idim) = inidim - 1
                    me%ibmx(idim) = inidim + 1

                    !  Distinguish the boundaries.
                    if (inidim==0) me%ibmn(idim) = 0
                    if (inidim==inmx(idim)) me%ibmx(idim) = inmx(idim)

                    !  Initialize the basis indices.
                    me%ib(idim) = me%ibmn(idim)
                end do

                !  Multiply by the extrapolation parameter (this acts as a
                !  smoothing weight).
                dcwght = swght*dcwght

                !  The COEF array serves as a row of the least squares matrix.
                !  Its value is zero except for columns corresponding to functions
                !  which are non-zero at the node.
                do icol = 1,ncol
                    coef(icol) = 0.0_wp
                end do

                !  The 2nd derivative of a function of MDIM variables may be thought
                !  of as a symmetric MDIM x MDIM matrix of 2nd order partial
                !  derivatives.  Traverse the upper triangle of this matrix and,
                !  for each element, compute a row of the least squares matrix.

                do idm = 1,me%mdim
                    do jdm = idm,me%mdim
                        do idim = 1,me%mdim
                            nderiv(idim) = 0
                        end do

                        boundary = .true.
                        !  Off-diagonal elements appear twice by symmetry, so the corresponding
                        !  row is weighted by a factor of 2.
                        rowwt = 2.0_wp*dcwght
                        if (jdm==idm) then
                            !  Diagonal.
                            rowwt = dcwght
                            nderiv(jdm) = 2
                            if (in(idm)/=0 .and. in(idm)/=inmx(idm)) then
                                boundary = .false.
                            end if
                        end if
                        if (boundary) then
                            !  Node is at boundary.
                            !
                            !  Normal 2nd derivative constraint at boundary is not appropriate for
                            !  natural splines (2nd derivative 0 by definition).  Substitute
                            !  a 1st derivative constraint.
                            nderiv(idm) = 1
                            nderiv(jdm) = 1
                        end if
                        irow = irow + 1

                        basis : do
                            !  Begining of basis index loop - traverse all indices corresponding
                            !  to basis functions which are non-zero at X.
                            !  The indices are in IB and are passed through common to BASCMP.
                            call me%bascmp(x,nderiv,xmin,nodes,icol,basm)

                            !  BASCMP computes ICOL and BASM where BASM is the value at X of the
                            !  N-dimensional basis function corresponding to column ICOL.
                            coef(icol) = rowwt*basm

                            !  Increment the basis indices.
                            do idim = 1,me%mdim
                                me%ib(idim) = me%ib(idim) + 1
                                if (me%ib(idim)<=me%ibmx(idim)) cycle basis
                                me%ib(idim) = me%ibmn(idim)
                            end do

                            !  End of basis index loop.
                            exit basis
                        end do basis

                        !  Send row of least squares matrix to reduction routine.
                        call me%suprls(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef,reserr,lserr)
                        if (lserr/=0) then
                            ierror = 107
                            call cfaerr(ierror, &
                                ' splcc or splcw - suprls failure '//&
                                '(this usually indicates insufficient input data)')
                        end if
                    end do
                end do

            end if

            !  Increment node indices.
            do idim = 1,me%mdim
                in(idim) = in(idim) + 1
                if (in(idim)<=inmx(idim)) cycle node_index
                in(idim) = 0
            end do

            exit node_index !  End of node index loop.

        end do node_index

    end if

    !  Call for least squares solution in COEF array.
    irow = 0
    call me%suprls(irow,coef,ncol,rhs,work(nwrk1),nwlft,coef,reserr,lserr)
    if (lserr/=0) then
        ierror = 107
        call cfaerr(ierror, &
            ' splcc or splcw - suprls failure '//&
            '(this usually indicates insufficient input data)')
    end if

end subroutine splcw
!*****************************************************************************************

!*****************************************************************************************
!>
!  N-dimensional cubic spline derivative evaluation.
!
!  A grid of evenly spaced nodes in `ndim` space is
!  defined by the arguments `xmin`, `xmax` and `nodes`.
!  A linear basis for the class of natural splines
!  on these nodes is formed, and to each basis
!  function corresponds a coefficient in the array
!  `coef` (computed in [[splcc]] or [[splcw]]).  Using
!  `nderiv` to indicate the appropriate partial
!  derivatives, each basis function is evaluated
!  at the point `x` in `ndim` space.  These values are
!  then multiplied by the corresponding
!  coefficient and summed to form the function
!  result.
!
!### See also
!  * [[splfe]]
!
!@note The original version of this routine would stop for an error.
!      Now it just returns.
!
!@note `coef`, `xmin`, `xmax` and `nodes` must be exactly
!      retained from the call to [[splcc]] (or [[splcw]]).

function splde(me,ndim,x,nderiv,coef,xmin,xmax,nodes,ierror)

    class(splpak_type),intent(inout) :: me
    real(wp) :: splde !! the function value returned is the partial
                      !! derivative (indicated by `nderiv`) of the
                      !! spline evaluated at `x`.
    integer,intent(in) :: ndim !! the dimensionality of the problem.  the
                               !! spline is a function of `ndim` variables or
                               !! coordinates and thus a point in the
                               !! independent variable space is an `ndim` vector.
                               !! `ndim` must be in the range `1 <= ndim <= 4`.
    real(wp),intent(in) :: x(ndim) !! an `ndim` vector describing the point in the
                                   !! independent variable space at which the
                                   !! spline is to be evaluated.
    real(wp),intent(out) :: coef(*) !! the array of coefficients which determine the
                                    !! spline.  each coefficient corresponds to a
                                    !! particular basis function which in turn
                                    !! corresponds to a node in the node grid.  this
                                    !! correspondence between the node grid and the
                                    !! array `coef` is as if `coef` were an
                                    !! ndim-dimensional fortran array with
                                    !! dimensions `nodes(1),...,nodes(ndim)`, i.e., to
                                    !! store the array linearly, the leftmost
                                    !! indices are incremented most frequently.
                                    !! coef may be computed by using routines [[splcc]]
                                    !! or [[splcw]].
                                    !!
                                    !! the dimension is assumed to be
                                    !! `coef(nodes(1)*...*nodes(ndim))`.
    real(wp),intent(in) :: xmin(ndim) !! a vector describing the lower extreme corner
                                      !! of the node grid.  a set of evenly spaced
                                      !! nodes is formed along each coordinate axis
                                      !! and `xmin(idim)` is the location of the first
                                      !! node along the `idim` axis.
    real(wp),intent(in) :: xmax(ndim) !! a vector describing the upper extreme corner
                                      !! of the node grid.  a set of evenly spaced
                                      !! nodes is formed along each coordinate axis
                                      !! and `xmax(idim)` is the location of the last
                                      !! node along the `idim` axis.
    integer,intent(in) :: nderiv(ndim) !! an `ndim` vector of integers specifying the
                                       !! partial derivative to be evaluated.  the
                                       !! order of the derivative along the `idim` axis
                                       !! is `nderiv(idim)`.  these integers must be in
                                       !! the range `0 <= nderiv(idim) <= 2`.
    integer,intent(in) :: nodes(ndim) !! a vector of integers describing the number of
                                      !! nodes along each axis.  `nodes(idim)` is the
                                      !! the number of nodes (counting endpoints)
                                      !! along the `idim` axis and determines the
                                      !! flexibility of the spline in that coordinate
                                      !! direction.  `nodes(idim)` must be >= 4 but
                                      !! may be as large as the arrays `coef` and `work`
                                      !! allow.
                                      !!
                                      !! *note:*  the node grid is completely defined by
                                      !! the arguments `xmin`, `xmax` and `nodes`.
                                      !! the spacing of this grid in the `idim`
                                      !! coordinate direction is
                                      !! `dx(idim) = (xmax(idim)-xmin(idim)) / (nodes(idim)-1)`.
                                      !! a node in this grid may be indexed by
                                      !! an `ndim` vector of integers
                                      !! `(in(1),...,in(ndim))` where
                                      !! `1 <= in(idim) <= nodes(idim)`.  the
                                      !! location of such a node may be
                                      !! represented by an `ndim` vector
                                      !! `(x(1),...,x(ndim)) ` where
                                      !! `x(idim) = xmin(idim)+(in(idim)-1) * dx(idim)`.
    integer,intent(out) :: ierror !! an error flag with the following meanings:
                                  !!
                                  !! *   0  no error.
                                  !! * 101  ndim is < 1 or is > 4.
                                  !! * 102  nodes(idim) is < 4 for some idim.
                                  !! * 103  xmin(idim) = xmax(idim) for some idim.
                                  !! * 104  nderiv(idim) is < 0 or is > 2 for some idim.

    real(wp) :: xrng,sum,basm
    integer :: iibmx,idim,nod,it,iib,icof

    ierror = 0
    me%mdim = ndim
    if (me%mdim<1) then
        ierror = 101
        call cfaerr(ierror, &
            ' splfe or splde - NDIM is less than 1')
        return
    end if
    iibmx = 1
    do idim = 1,me%mdim
        nod = nodes(idim)
        if (nod<4) then
            ierror = 102
            call cfaerr(ierror, &
                ' splfe or splde - NODES(IDIM) is less than  4for some IDIM')
            return
        end if
        xrng = xmax(idim) - xmin(idim)
        if (xrng==0.0_wp) then
            ierror = 103
            call cfaerr(ierror, &
                ' splfe or splde - XMIN(IDIM) = XMAX(IDIM) for some IDIM')
            return
        end if
        if (nderiv(idim)<0 .or. nderiv(idim)>2) then
            ierror = 104
            call cfaerr(ierror, &
                ' splde - NDERIV(IDIM) IS less than 0 or greater than 2 for some IDIM')
        end if

        !  DX(IDIM) is the node spacing along the IDIM coordinate.
        me%dx(idim) = xrng/real(nod-1,wp)
        me%dxin(idim) = 1.0_wp/me%dx(idim)

        !  Compute indices of basis functions which are nonzero at X.
        it = me%dxin(idim)*(x(idim)-xmin(idim))

        !  IBMN must be in the range 0 to NODES-2.
        me%ibmn(idim) = min(max(it-1,0),nod-2)

        !  IBMX must be in the range 1 to NODES-1.
        me%ibmx(idim) = max(min(it+2,nod-1),1)
        iibmx = iibmx* (me%ibmx(idim)-me%ibmn(idim)+1)
        me%ib(idim) = me%ibmn(idim)
    end do

    sum = 0.0_wp
    iib = 0

    basis_index : do
        !  Begining of basis index loop - traverse all indices corresponding
        !  to basis functions which are nonzero at X.
        iib = iib + 1

        !  The indices are in IB and are passed through common to BASCMP.
        call me%bascmp(x,nderiv,xmin,nodes,icof,basm)

        !  BASCMP computes ICOF and BASM where BASM is the value at X of the
        !  N-dimensional basis function corresponding to COEF(ICOF).
        sum = sum + coef(icof)*basm
        if (iib<iibmx) then
            !  Increment the basis indices.
            do idim = 1,me%mdim
                me%ib(idim) = me%ib(idim) + 1
                if (me%ib(idim)<=me%ibmx(idim)) cycle basis_index
                me%ib(idim) = me%ibmn(idim)
            end do
        end if

        exit basis_index !  End of basis index loop.
    end do basis_index

    splde = sum

end function splde
!*****************************************************************************************

!*****************************************************************************************
!>
!  N-dimensional cubic spline function evaluation.
!
!  Except for lack of derivative capability, this
!  function is identical to function [[splde]] in
!  usage.  The argument list is also identical
!  except for the omission of `nderiv`.
!
!### See also
!  * [[splde]]
!
!@note `coef`, `xmin`, `xmax` and `nodes` must be exactly
!      retained from the call to [[splcc]] (or [[splcw]]).

function splfe(me,ndim,x,coef,xmin,xmax,nodes,ierror)

    class(splpak_type),intent(inout) :: me
    real(wp) :: splfe
    integer,intent(in) :: ndim
    real(wp),intent(in) :: x(ndim)
    real(wp),intent(out) :: coef(*)
    real(wp),intent(in) :: xmin(ndim)
    real(wp),intent(in) :: xmax(ndim)
    integer,intent(in) :: nodes(ndim)
    integer,intent(out) :: ierror

    integer,dimension(ndim) :: nderiv

    nderiv = 0
    splfe = me%splde(ndim,x,nderiv,coef,xmin,xmax,nodes,ierror)

end function splfe
!*****************************************************************************************

!*****************************************************************************************
!>
!  To determine the least squares solution of a
!  large overdetermined linear system.
!
!  Given the `m` by `n` matrix `r` (`m >= n`)
!  and the `m`-vector `b`,
!  this routine calculates the `n`-vector `x` such
!  that the euclidean norm of the residue (`r*x-b`)
!  is minimized.  the subroutine accepts rows of
!  the matrix one by one so that the entire matrix
!  need not be stored at one time.  this allows
!  large problems to be solved without peripheral
!  storage.  the length of the rows is limited by
!  the amount of scratch storage which can be set
!  aside for use by the routine.  there is no
!  restriction on the number of rows.
!
!### Usage
!  [[suprls]] is called once for
!  each row of the matrix.  A final call returns
!  the solution vector and the euclidean norm
!  of the residual. This following sequence would
!  process the `m` by `n` matrix `r` and the right hand
!  side `m`-vector `b`
!
!```fortran
!  do i = 1,m
!    do j = 1,n
!      ! here set rowi(j) to the (i,j) element of r
!    end do
!    ! here set bi to the ith component of b.
!    call suprls(i,rowi,n,bi,a,nn,soln,err,ier)
!  end do
!  call suprls (0,rowi,n,bi,a,nn,soln,err,ier)
!```
!
!### Algorithm
!  given the `m` by `n` matrix `r` (`m>=n`) and the
!  `m`-vector `b`, we wish to find an `n`-vector `x`
!  such that
!```
!      e = l2-norm of  r*x-b
!```
!  is minimized.  since the euclidean norm is
!  invariant under orthogonal transformation,
!  `r` and `b` may be premultiplied by any
!  orthogonal matrix without changing the
!  norm of the residual (`r*x-b`).  `r` is reduced
!  to upper triangular form by premultiplying
!  `r` and `b` by a sequence of householder and
!  rotation matrices.  when the reduction is
!  complete, the norm of the residual takes the
!  form
!```
!    e =  l2 norm(t*x-b(n))+l2 norm(b(m-n))
!```
!  where `t` is an `n` by `n` upper triangular
!  matrix, `b(n)` is a vector of the first `n`
!  components of `b`, `b(m-n)` is a vector of
!  the remaining `(m-n)` components of `b`.  `e` is
!  minimized by taking `x` to be the solution
!  of the system  `t*x=b(n)`.  this triangular
!  system is therefore solved to give the
!  required least squares solution.  the norm
!  of the residual is then the l2-norm of `b(m-n)`.
!
!  at each phase of the reduction, as many rows
!  as space permits are entered into the scratch
!  area.  householder transformations are then
!  used to zero out subdiagonal elements.  space
!  is saved by eliminating storage for the
!  zero subdiagonal terms.  if there is room
!  for only one new row, rotation rather than
!  householder matrices are used for greater
!  speed.  when all `m` rows have been entered,
!  reduction is completed and the triangular
!  system solved.
!
!### Reference
!  * Hanson, R.J., and Lawson, C.L., "Extensions and
!  applications of the householder algorithm
!  for solving linear least squares problems".
!  Math. of comp. vol.23, pp. 787-812. (1969)
!
!### Accuracy
!  This will depend upon the size and condition
!  of the matrix.  Near machine accuracy may be
!  expected for well conditioned systems of
!  moderate size.  If ill conditioning is
!  suspect, a version using pivoting may be
!  necessary.
!
!### History
!  * Original FORTRAN 66 version written in May 1972 by
!    A.K. Cline of NCAR'S Scientific Computing Division.

subroutine suprls(me,i,rowi,n,bi,a,nn,soln,err,ier)

    class(splpak_type),intent(inout) :: me
    integer,intent(in) :: i !! the index of the row being entered.  (`i` is 1
                            !! for the first call, increases by 1 for each
                            !! call, and is `m` when the final row is
                            !! entered).  After the final row has been
                            !! entered, [[suprls]] is called with `i = 0` to
                            !! complete the reduction and solution.
    integer,intent(in) :: nn !! length of scratch array `a`. `nn` must be at
                             !! least `n*(n+5)/2+1`. For speed, `nn` should be
                             !! as large as possible up to a maximum of
                             !! `(n+1)*m`.
    integer,intent(in) :: n !! the length of the rows of the matrix (i.e.,
                            !! the number of columns). `n <= m`, where `m` is
                            !! the number of rows.
    real(wp),intent(in) :: rowi(n) !! a vector which on the `i`th call contains the `n`
                                   !! components of the `i`th row of the matrix.  the
                                   !! dimension of `rowi` in calling program must be
                                   !! at least `n`.
    real(wp),intent(in) :: bi !! on the `i`th call, `bi` contains the `i`th element
                              !! of the right hand side vector `b`.
    real(wp),intent(inout) :: a(nn) !! a working array which must not be changed
                                    !! between the successive calls to [[suprls]].
    real(wp),intent(out) :: soln(n) !! the `n`-components of the solution vector are
                                    !! returned in this array after the final call
                                    !! to [[suprls]].
    real(wp),intent(out) :: err !! the euclidean norm of the residual is
                                !! returned in `err` after the final call to
                                !! [[suprls]].
    integer,intent(out) :: ier !! error parameter.
                               !! fatal errors:
                               !!
                               !!  * 32 -- insufficient scratch storage provided,
                               !!    must have `nn >= n*(n+5)/2+1`.
                               !!  * 33 -- array has too few rows.  must have
                               !!    `m >= n`.
                               !!  * 34 -- system is singular.
                               !!  * 35 -- values of `i` not in sequence.

    real(wp) :: s,temp,temp1,cn,sn
    integer :: j,ilj,ilnp,nreq,k,&
               idiag,i1,i2,ii,jp1,lmkm1,j1,jdel,idj,iijd,&
               i1jd,k11,k1m1,i11,np1mk,lmk,imov,iii,iiim,iim1,&
               ilk,npk,ilii,npii
    logical :: complete_reduction !! Routine entered with `I<=0` means complete
                                  !! the reduction and store the solution in `SOLN`.

    real(wp),parameter :: tol = 1.0e-18_wp !! small number tolerance

    ier = 0
    complete_reduction = i <= 0

    if (.not. complete_reduction) then

        if (i<=1) then

            !  Set up quantities on first call.
            me%iold = 0
            me%np1 = n + 1

            !  Compute how many rows can be input now.
            me%l = nn/me%np1
            me%ilast = 0
            me%il1 = 0
            me%k = 0
            me%k1 = 0
            me%errsum = 0.0_wp
            nreq = ((n+5)*n+2)/2

            !  Error exit if insufficient scratch storage provided.
            if (nn<nreq) then
                ier = 32
                write(*,*) 'nn   = ', nn
                write(*,*) 'nreq = ', nreq
                call cfaerr(ier, &
                            ' suprls - insufficient scratch storage provided. '//&
                            'at least ((N+5)*N+2)/2 locations needed')
                return
            end if

        end if

        !  Error exit if (I-IOLD)/=1.
        if ((i-me%iold)/=1) then
            ier = 35
            write(*,*) 'i    =',i
            write(*,*) 'me%iold =',me%iold
            call cfaerr(ier,' suprls - values of I not in sequence')
            return
        end if

        !  Store the row in the scratch storage.
        me%iold = i
        do j = 1,n
            ilj = me%ilast + j
            a(ilj) = rowi(j)
        end do
        ilnp = me%ilast + me%np1
        a(ilnp) = bi
        me%ilast = me%ilast + me%np1
        me%isav = i
        if (i<me%l) return

    end if

    main : do

        if (.not. complete_reduction) then

            if (me%k/=0) then
                me%k1 = min(me%k,n)
                idiag = -me%np1
                if (me%l-me%k==1) then
                    !  Apply rotations to zero out the single new row.
                    do j = 1,me%k1
                        idiag = idiag + (me%np1-j+2)
                        i1 = me%il1 + j
                        if (abs(a(i1))<=tol) then
                            s = sqrt(a(idiag)*a(idiag))
                        else if (abs(a(idiag))<tol) then
                            s = sqrt(a(i1)*a(i1))
                        else
                            s = sqrt(a(idiag)*a(idiag)+a(i1)*a(i1))
                        end if
                        if (s==0.0_wp) cycle
                        temp = a(idiag)
                        a(idiag) = s
                        s = 1.0_wp/s
                        cn = temp*s
                        sn = a(i1)*s
                        jp1 = j + 1
                        do j1 = jp1,me%np1
                            jdel = j1 - j
                            idj = idiag + jdel
                            temp = a(idj)
                            i1jd = i1 + jdel
                            a(idj) = cn*temp + sn*a(i1jd)
                            a(i1jd) = -sn*temp + cn*a(i1jd)
                        end do
                    end do
                else
                    !  Apply householder transformations to zero out new rows.
                    do j = 1,me%k1
                        idiag = idiag + (me%np1-j+2)
                        i1 = me%il1 + j
                        i2 = i1 + me%np1* (me%l-me%k-1)
                        s = a(idiag)*a(idiag)
                        do ii = i1,i2,me%np1
                            s = s + a(ii)*a(ii)
                        end do
                        if (s==0.0_wp) cycle
                        temp = a(idiag)
                        a(idiag) = sqrt(s)
                        if (temp>0.0_wp) a(idiag) = -a(idiag)
                        temp = temp - a(idiag)
                        temp1 = 1.0_wp / (temp*a(idiag))
                        jp1 = j + 1
                        do j1 = jp1,me%np1
                            jdel = j1 - j
                            idj = idiag + jdel
                            s = temp*a(idj)
                            do ii = i1,i2,me%np1
                                iijd = ii + jdel
                                s = s + a(ii)*a(iijd)
                            end do
                            s = s*temp1
                            a(idj) = a(idj) + s*temp
                            do ii = i1,i2,me%np1
                                iijd = ii + jdel
                                a(iijd) = a(iijd) + s*a(ii)
                            end do
                        end do
                    end do
                end if

                if (me%k>=n) then
                    lmkm1 = me%l - me%k

                    !  Accumulate residual sum of squares.
                    do ii = 1,lmkm1
                        ilnp = me%il1 + ii*me%np1
                        me%errsum = me%errsum + a(ilnp)*a(ilnp)
                    end do
                    if (i<=0) exit main
                    me%k = me%l
                    me%ilast = me%il1

                    !  Determine how many new rows may be input on next iteration.
                    me%l = me%k + (nn-me%ilast)/me%np1
                    return
                end if
            end if

            k11 = me%k1 + 1
            me%k1 = min(me%l,n)
            if (me%l-me%k/=1) then
                k1m1 = me%k1 - 1
                if (me%l>n) k1m1 = n
                i1 = me%il1 + k11 - me%np1 - 1

                !  Perform householder transformations to reduce rows to upper
                !  triangular form.
                do j = k11,k1m1
                    i1 = i1 + (me%np1+1)
                    i2 = i1 + (me%l-j)*me%np1
                    s = 0.0_wp
                    do ii = i1,i2,me%np1
                        s = s + a(ii)*a(ii)
                    end do
                    if (s==0.0_wp) cycle
                    temp = a(i1)
                    a(i1) = sqrt(s)
                    if (temp>0.0_wp) a(i1) = -a(i1)
                    temp = temp - a(i1)
                    temp1 = 1.0_wp/ (temp*a(i1))
                    jp1 = j + 1
                    i11 = i1 + me%np1
                    do j1 = jp1,me%np1
                        jdel = j1 - j
                        i1jd = i1 + jdel
                        s = temp*a(i1jd)
                        do ii = i11,i2,me%np1
                            iijd = ii + jdel
                            s = s + a(ii)*a(iijd)
                        end do
                        s = s*temp1
                        i1jd = i1 + jdel
                        a(i1jd) = a(i1jd) + s*temp
                        do ii = i11,i2,me%np1
                            iijd = ii + jdel
                            a(iijd) = a(iijd) + s*a(ii)
                        end do
                    end do
                end do
                if (me%l>n) then
                    np1mk = me%np1 - me%k
                    lmk = me%l - me%k
                    ! Accumulate residual sum of squares.
                    do ii = np1mk,lmk
                        ilnp = me%il1 + ii*me%np1
                        me%errsum = me%errsum + a(ilnp)*a(ilnp)
                    end do
                end if
            end if
            imov = 0
            i1 = me%il1 + k11 - me%np1 - 1

            !  Squeeze the unnecessary elements out of scratch storage to
            !  allow space for more rows.
            do ii = k11,me%k1
                imov = imov + (ii-1)
                i1 = i1 + me%np1 + 1
                i2 = i1 + me%np1 - ii
                do iii = i1,i2
                    iiim = iii - imov
                    a(iiim) = a(iii)
                end do
            end do
            me%ilast = i2 - imov
            me%il1 = me%ilast
            if (i<=0) exit main
            me%k = me%l

            !  Determine how many new rows may be input on next iteration.
            me%l = me%k + (nn-me%ilast)/me%np1
            return

        end if

        ! Complete reduction and store solution in SOLN.
        complete_reduction = .false. ! reset this flag
        me%l = me%isav

        ! Error exit if L less than N.
        if (me%l<n) then
            ier = 33
            call cfaerr(ier,' suprls - array has too few rows.')
            return
        end if

        ! K/=ISAV means further reduction needed.
        if (me%k==me%isav) exit main

    end do main

    me%ilast = (me%np1* (me%np1+1))/2 - 1
    if (a(me%ilast-1)==0.0_wp) then
        ! Error return if system is singular.
        ier = 34
        call cfaerr(ier,' suprls - system is singular.')
        return
    end if

    ! Solve triangular system into ROWI.
    soln(n) = a(me%ilast)/a(me%ilast-1)
    do ii = 2,n
        iim1 = ii - 1
        me%ilast = me%ilast - ii
        s = a(me%ilast)
        do k = 1,iim1
            ilk = me%ilast - k
            npk = me%np1 - k
            s = s - a(ilk)*soln(npk)
        end do
        me%k = k ! JW : is this necessary ???
        ilii = me%ilast - ii
        if (a(ilii)==0.0_wp) then
            ! Error return if system is singular.
            ier = 34
            call cfaerr(ier,' suprls - system is singular.')
            return
        end if
        npii = me%np1 - ii
        soln(npii) = s/a(ilii)
    end do

    ! Store residual norm.
    err = sqrt(me%errsum)

end subroutine suprls
!*****************************************************************************************

!*****************************************************************************************
    end module splpak_module
!*****************************************************************************************