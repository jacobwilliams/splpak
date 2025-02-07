!*****************************************************************************************
!>
!  Units test for splpak.
!  Test a linear function

    program splpak_test_linesr

    use splpak_module, wp => splpak_wp
    use iso_fortran_env
    use pyplot_module

    implicit none

    integer,parameter :: ndim = 1    !! 1d problem
    integer,parameter :: nxdata = 20 !! number of points in x
    integer,dimension(ndim),parameter :: nodes = [10]
    integer,parameter :: ncol = product(nodes)
    integer,parameter :: nwrk = ncol*(ncol+1) + 1   ! is doc wrong? do we need the + 1 (otherwise, it fails)...
    integer,parameter :: ncf = ncol
    integer,parameter :: nxdata_est = 100 !! number of points for estimate plot

    real(wp),dimension(ndim,nxdata) :: xdata
    real(wp),dimension(nxdata) :: ydata
    real(wp),dimension(ndim)   :: xmin
    real(wp),dimension(ndim)   :: xmax
    real(wp),dimension(nwrk)   :: work
    real(wp),dimension(ncf)    :: coef
    real(wp),dimension(nxdata) :: wdata !! weights
    real(wp),dimension(nxdata_est) :: xdata_est, ydata_est
    real(wp),dimension(ndim)   :: x

    integer :: ierror, istat
    integer :: i !! counter
    real(wp) :: xtrap
    real(wp) :: tru, err, errmax, f, fleft, fright
    type(pyplot) :: plt
    character(len=10) :: nodes_str !! string version of `nodes`
    type(splpak_type) :: solver, solver2
    integer,dimension(2),parameter :: figsize = [20,10] !! figure size for plot

    xtrap = 1.0_wp
    xmin(1) = 0.0_wp
    xmax(1) = 1.0_wp
    do i=1,nxdata
        wdata(i) = 1.0_wp   ! weight
        xdata(1,i) = real(i-1,wp)/real(nxdata-1,wp)
        ydata(i) = f1(xdata(1,i))
    end do

    ! write(*,'(a,*(f8.3,","))') 'xdata = ', xdata
    ! write(*,'(a,*(f8.3,","))') 'ydata = ', ydata

    ! initialize:
    call solver%initialize(1,xdata,1,ydata,wdata,nxdata,xmin,xmax, &
                           nodes,xtrap,coef,ncf,work,nwrk,ierror)

    write(*,*) 'splcw ierror = ', ierror
    if (ierror /= 0) error stop 'error calling splcw'

    ! compute max error at interpolation points
    errmax = 0.0_wp
    do i=1,nxdata_est
        xdata_est(i) = real(i-1,wp)/nxdata_est
        x(1) = xdata_est(i)
        f = solver%evaluate(ndim,x,coef,xmin,xmax,nodes,ierror)
        ydata_est(i) = f
        if (ierror /= 0) error stop 'error calling splfe'
        tru    = f1(xdata_est(i))
        err    = abs(tru-f)
        errmax = max(err,errmax)
    end do
    write(*,*) 'splfe errmax [linear] = ', errmax
    if (abs(errmax)>1.0e-1_wp) error stop 'errmax too large'

    ! write(*,'(a,*(f8.3,","))') 'xdata_est = ', xdata_est
    ! write(*,'(a,*(f8.3,","))') 'ydata_est = ', ydata_est

    ! test the derivatives also:
    fleft = solver%evaluate(ndim,[0.0_wp],[1],coef,xmin,xmax,nodes,ierror)
    if (ierror /= 0) error stop 'error calling splde'
    errmax = fleft - 2.0_wp
    write(*,*) 'splde errmax left [linear] = ', errmax
    if (abs(errmax)>1.0e-12_wp) error stop 'errmax too large'

    fright = solver%evaluate(ndim,[1.0_wp],[1],coef,xmin,xmax,nodes,ierror)
    if (ierror /= 0) error stop 'error calling splde'
    errmax = fleft - 2.0_wp
    write(*,*) 'splde errmax right [linear] = ', errmax
    if (abs(errmax)>1.0e-12_wp) error stop 'errmax too large'

    write(nodes_str,'(I10)') nodes(1); nodes_str = adjustl(nodes_str)

    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                        figsize=figsize,font_size=20,axes_labelsize=20,&
                        xtick_labelsize=20, ytick_labelsize=20,&
                        legend_fontsize=20,&
                        title='splpak_test_linear',legend=.true.)
    call plt%add_plot(xdata(1,:),ydata,&
                        label='Original points',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(xdata_est,ydata_est,&
                        label='Least squares bspline with '//trim(nodes_str)//' nodes',&
                        linestyle='r-',markersize=2,linewidth=2,istat=istat)
    call plt%savefig(pyfile='splpak_test.py', figfile='splpak_test_linear.png',istat=istat)

    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        f1 = 2.0_wp * x
        end function f1

    end program splpak_test_linesr