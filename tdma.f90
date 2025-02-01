 program compact_fd 
    implicit none
    
    integer , parameter :: dp = selected_real_kind(15) !
    integer :: i                        !loop indices
    integer , parameter :: n = 51   ! no of grid points
    real(dp) :: h                        !grid spacing 
    real(dp) , parameter :: pi = 4.0_dp*atan(1.0_dp)
    real(dp) , dimension(n) :: x , f , df_exact , df_num , error , error_percent! arrays for storing grid , function, derivatives
    
    h = 2.0_dp*pi/(n-1) ! grid spacing 
    
    !grid and function initialization
    do i = 1, n 
        x(i) = (i-1)*h 
        f(i) = sin(x(i))
        df_exact(i) = cos(x(i)) !exact derivative
    end do 
    
    call compact_fd6(n, f, df_num , h )
    
    error = abs(df_num - df_exact)
    
    do i = 1, n
   if (abs(df_exact(i)) < tiny(1.0_dp)) then
    error_percent(i) = 0.0_dp  ! Avoid division by zero
  else
    error_percent(i) = 100.0_dp * error(i) / abs(df_exact(i))
  end if
end do
    
    
    open(unit=10 , file='results.txt' , status='replace')
    write(10 , '(A)') 'x, f(x) , df_exact(x) , df_num(x) , error '
    do i=1,n 
        write(10 , '(6E16.8)') x(i) , f(i) , df_exact(i) , df_num(i) , error(i) 
    end do 
    
    close(10)
    
    print * , 'results written to results.txt'
    print *, 'maximum error:' , maxval(error)
    !print * , 'maximum error percentage:' , maxval(error_percent)
    
    
    contains 
    
    
    subroutine compact_fd6(n , f , df_num  , h )
    
    implicit none 
    integer , intent(in) :: n 
    integer :: i ,nt
    real(dp) :: alpha , beta , gamma
    real(dp) , allocatable :: a(:) , b(:) , c(:) , d(:) 
    real(dp) , intent(in) :: h
    real(dp) , intent(out) , dimension(n) :: df_num
    real(dp) , dimension(n-4) :: df_num_sub
    real(dp) , intent(in) , dimension(n) :: f 
    
    nt = n-4
    !coefficients for interior points 
    alpha = 1.0_dp/3.0_dp
    beta = 14.0_dp/18.0_dp
    gamma = 1.0_dp/36.0_dp
    
 ! explicitly calculating using one-sided 4th order accurate schemes   
    
! For i=1 (4th order forward)
df_num(1) = (-25.0_dp*f(1) + 48.0_dp*f(2) - 36.0_dp*f(3) + 16.0_dp*f(4) - 3.0_dp*f(5)) / (12.0_dp*h)
! For i=2 (4th order)
df_num(2) = (-3.0_dp*f(1) - 10.0_dp*f(2) + 18.0_dp*f(3) - 6.0_dp*f(4) + 1.0_dp*f(5)) / (12.0_dp*h)
! For i=n-1 (4th order)
df_num(n-1) = (3.0_dp*f(n) + 10.0_dp*f(n-1) - 18.0_dp*f(n-2) + 6.0_dp*f(n-3) - 1.0_dp*f(n-4)) / (12.0_dp*h)
! For i=n (4th order backward)
df_num(n) = (25.0_dp*f(n) - 48.0_dp*f(n-1) + 36.0_dp*f(n-2) - 16.0_dp*f(n-3) + 3.0_dp*f(n-4)) / (12.0_dp*h)
   
! allocate arrays for tridiagonal matrix(i=3,n-2) and rhs
    allocate(a(nt-1) , b(nt) , c(nt-1) , d(nt) ) 
    
    ! initialize tri diagonal matrix and rhs(d)
    
    a(:) = alpha
    b(:) = 1.0_dp
    c(:) = alpha
    
    
    do i=1,nt
        d(i) = (beta*(f(i+3) - f(i+1)) + gamma*(f(i+4) - f(i)))/h
    end do 

call tdma(a, b, c, d, df_num_sub, nt)
!copying result back to original array 
df_num(3:n-2) = df_num_sub

deallocate(a , b , c , d)
end subroutine compact_fd6

subroutine tdma(a , b, c ,d , x , nt)
    
implicit none

integer :: i
integer , intent(in) :: nt 
real(dp) , intent(in) :: a(nt-1) , b(nt) ,c(nt-1) , d(nt)
real(dp) , intent(out) :: x(nt)
real(dp) , dimension(nt) :: e, f
real(dp) :: t

   
    
    !forward sweep
    e(1) = c(1)/b(1)
    f(1) = d(1)/b(1)
    
    
    do i = 2, nt-1
        t = b(i) - a(i-1)*e(i-1)
        e(i) = c(i)/t
        f(i) = (d(i)-a(i)*f(i-1))/t
    end do 
    
    !explicitly calculating for i=nt(to skip unnecesaary calculations for e(nt))
    t = b(nt) - a(nt-1)*e(nt-1)
    f(nt) = (d(nt) - a(nt-1)*f(nt-1))/t
    
    !backward sweep 
    
    x(nt) = f(nt) 
    do i = nt-1 , 1 , -1
        x(i) = f(i) - e(i)*x(i+1)
    end do 
    end subroutine tdma

    end program compact_fd
    
