program test
    implicit none

    real*8, dimension(3,5) ::  weno_coeffs_de_avg
    integer :: i, j
    weno_coeffs_de_avg = reshape((/ -1.0/24.0, 0.,0.,&
                                    1.0/12.0, -1.0/24.0, 0.,&
                                    23.0/24.0, 13.0/12.0, 23.0/24.0,&
                                    0.,-1.0/24.0,1.0/12.0,&
                                    0.,0.,-1.0/24.0&
    /), shape(weno_coeffs_de_avg))

    loopy: do j = 1, 3
        loopx: do i = 1, 5
            write(*,*) "x=", i, ",y=", j, "=?=", weno_coeffs_de_avg(j, i)
        end do loopx
    end do loopy
    
end program test