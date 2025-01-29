program compact_fd 
    implicit none
    
    integer , parameter :: dp = selected_real_kind(15) !
    integer :: i                        !loop indices
    integer , parameter :: n = 101      ! no of grid points
    real(dp) :: h                        !grid spacing 
    real(dp) , parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp) , dimension(n) :: x , f ,df_exact , df_ num , error ! arrays for storing grid , function, derivatives
    
    h = 2.0_dp*pi/(n-1) ! grid spacing 
    
    !grid and function initialization
    do i = 1, n 
        x(i) = (i-1)*h 
        f(i) = sin(x(i))
        df_exact(i) = cos(x(i)) !exact derivative
    end do 
    
    call compact_fd6(n, f, df_num , h )
    
    error = abs(df_num - df_exact)
    
    open(unit=10 , file='results.txt' , status='replace'
    write(10 , '(A)') 'x, f(x) , df_exact(x) , df_num(x) , error'
    do i=1,n 
        write(10 , '(5F16.8)') x(i) , f(i) , df_exact(i) , df_num(i) , error(i)
    end do 
    
    close(10)
    
    print * , 'results written to results.txt'
    print *, 'maximum error:' , maxval(error)
    
    
    contains 
    
    
    subroutine compact_fd6(n , f  df_num  , h )
    
    implicit none 
    integer , intent(in) :: n 
    real(dp) , intent(in) :: h
    real(dp) , intent(out) , dimension(n) :: df_num
    real(dp) , intent(in) , dimension(n) :: f 
    !coeffcicents 
    
    l = 1.0_dp/3.0_dp
    m = 14.0_dp/18.0_dp
    n = 1.0_dp/36.0_dp
    
    ! allocate arrays for tridiagonal matrix and rhs
    allocate(a(n-1) , b(n) , c(n-1) , d(n) ) 
    
    ! initialize tri diagonal matrix and rhs 
    
    a(:) = l
    b(:) = 1.0_dp
    c(:) = l 
    d(:) = 0.0_dp
    
    do i=3,n-2
        d(i) = (a*(f(i+1) - f(i-1)) + b*(f(i+2) - f(i-2)))/h
    end do 
    
    
    !the boundary conditions lead to some A or C coefficients being zero, but the tridiagonal structure remains, so the algorithm works.
    !This method is computationally efficient, requiring O(n) operations instead of the O(n³) operations of Gaussian elimination for a full matrix.    
    !The algorithm’s stability depends on avoiding division by zero, but boundary treatments in the code ensure this.
    !This flexibility makes the Thomas algorithm ideal for compact finite difference schemes with mixed boundary and interior stencils.
!boundary points (one sided schemes)
d(1) = (-11.0_dp*f(1) + 18.0_dp*f(2) -9.0_dp*f(3) + 2.0_dp*f(4))/(6.0_dp*h)
d(2) = (-2.0_dp*f(1) - 3.0_dp*f(2) + 6.0_dp*f(3) - 1.0_dp*f(4))/(6.0_dp*h)
d(n-1) = (2.0_dp*f(n) + 3.0_dp*f(n-1) -6.0_dp*f(n-2) + 1.0_dp*f(n-2))/(6.0_dp*h)
d(n) = (11.0_dp*f(n) - 18.0_dp*f(n-1) + 9.0_dp*f(n-2) - 2.0_dp*f(n-3))/(6.0_dp*h)

call tdma(a , b  , c, d , df_num , n )

    subroutine tdma(a , b, c ,d , x , n)
    

    implicit none

    integer , parameter :: dp = selected_real_kind(15)
    integer :: i
    integer , intent(in) :: n 
    real(dp) , intent(in) : a(n-1) , b(n) ,c(n-1) , d(n)
    real(dp) , intent(out) : x
    real(dp) , dimension(n) : e, f
    real(dp) :: t
    
    !forward sweep
    e(1) = c(1)/b(1)
    f(1) = d(1)/b(1)
    
    do i = 2, n
        t = b(i) - a(i)*e(i-1)
        e(i) = c(i)/t
        f(i) = (d(i)-a(i)*f(i-1))/t
    end do 
    
    !backward sweep 
    
    x(n) = f(n) 
    do i = n-1 , 1 , -1
        x(i) = f(i) - e(i)*x(i+1)
    end do 
    end subroutine tdma

    end program compact_fd
    
