!*****************************************************************************************
!>
!  Units test for splpak.
!
!### See also
!  * `bspline_test` in [bspline-fortran](https://github.com/jacobwilliams/bspline-fortran)

    program splpak_test

    use splpak_module
    use iso_fortran_env, only: wp => real64, ip => int32
    use pyplot_module

    implicit none

    integer,parameter :: ndim = 1    !! 1d problem
    integer(ip),parameter :: nxdata = 20 !! number of points in x
    real(wp),dimension(1),parameter :: wdata = [-1.0_wp] !! don't use weights
    integer,dimension(ndim),parameter :: nodes = [10]
    integer,parameter :: ncol = product(nodes)
    integer,parameter :: nwrk = ncol*(ncol+1) + 1   ! is doc wrong? do we need the + 1 (otherwise, it fails)...
    integer,parameter :: ncf = ncol

    real(wp),dimension(ndim,nxdata) :: xdata
    real(wp),dimension(nxdata) :: ydata, ydata_est
    real(wp),dimension(ndim)   :: xmin
    real(wp),dimension(ndim)   :: xmax
    real(wp),dimension(nwrk)   :: work
    real(wp),dimension(ncf)    :: coef
    integer :: ierror, istat
    integer :: i !! counter
    real(wp) :: xtrap
    real(wp) :: tru, err, errmax, f
    type(pyplot) :: plt
    integer(ip),dimension(:),allocatable :: iseed
    real(wp) :: r !! random number
    integer(ip) :: isize !! for `random_seed`
    character(len=10) :: nodes_str !! string version of `nodes`

    call random_seed(size=isize)
    allocate(iseed(isize)); iseed = 42_ip
    call random_seed(put=iseed)

    xtrap = 1.0_wp
    xmin(1) = 0.0_wp
    xmax(1) = 1.0_wp
    do i=1,nxdata
        call random_number(r)
        r = (r - 0.5_wp)/ 10.0_wp   ! some random noise

        xdata(1,i) = real(i-1,wp)/real(nxdata-1,wp)
        ydata(i) = f1(xdata(1,i)) + r
    end do

    write(*,'(a,*(f8.3,","))') 'xdata = ', xdata
    write(*,'(a,*(f8.3,","))') 'ydata = ', ydata

    ! initialize:
    call splcw(1,xdata,1,ydata,wdata,nxdata,xmin,xmax, &
               nodes,xtrap,coef,ncf,work,nwrk,ierror)

    write(*,*) 'splcw ierror = ', ierror
    if (ierror /= 0) error stop 'error calling splcw'

    ! compute max error at interpolation points
    errmax = 0.0_wp
    do i=1,nxdata
        f = splfe (ndim,xdata(1,i),coef,xmin,xmax,nodes,ierror)
        ydata_est(i) = f
        if (ierror /= 0) error stop 'error calling splfe'
        tru    = f1(xdata(1,i))
        err    = abs(tru-f)
        errmax = max(err,errmax)
    end do
    write(*,*) 'splfe errmax = ', errmax
    if (abs(errmax)>1.0e-1_wp) error stop 'errmax too large'

    write(nodes_str,'(I10)') nodes(1); nodes_str = adjustl(nodes_str)

    call plt%initialize(grid=.true.,xlabel='x',ylabel='y',&
                        figsize=[20,10],font_size=20,axes_labelsize=20,&
                        xtick_labelsize=20, ytick_labelsize=20,&
                        legend_fontsize=20,&
                        title='splpak_test',legend=.true.)
    call plt%add_plot(xdata(1,:),ydata,&
                        label='Original points',&
                        linestyle='ko',markersize=5,linewidth=2,istat=istat)
    call plt%add_plot(xdata(1,:),ydata_est,&
                        label='Least squares bspline with '//trim(nodes_str)//' nodes',&
                        linestyle='r-',markersize=2,linewidth=2,istat=istat)
    call plt%savefig(pyfile='splpak_test.py', figfile='splpak_test.png',istat=istat)


    contains

        real(wp) function f1(x) !! 1d test function
        implicit none
        real(wp) :: x
        f1 = 0.5_wp * (x*exp(-x) + sin(x) )
        end function f1

    end program splpak_test