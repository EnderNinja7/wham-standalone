program weno
    implicit none
    !
    !                                               WHAM STANDALONE CODE
    !   PURPOSE: IMPLEMENT 3-DIMENSIONAL DE-AVERAGING STEPS OVER FORTRAN, BEFORE APPLYING THEM TO GRHYDRO_CON2PRIM
    !                                             AUTHOR: DHRUV SRIVASTAVA
    !
    !

    !Declare parameters and variables
    real*8, dimension(10,10,10)::data_
    integer :: i, j, k, itracer, nx, ny, nz, GRHydro_stencil
    real*8 x0, y0, z0, dxyz, x, y, z, weno_eps

    real*8, dimension(10, 10, 10):: temp1, temp2, a_center
    
    nx = 10
    ny = 10
    nz = 10
    x0 = 1.
    y0 = 1.
    z0 = 1.
    dxyz = 0.1
    GRHydro_stencil = 3
    weno_eps = 1e-10

    

    

    !Fill array with dummy data 

    loopx: do i = 1, nx
        loopy: do j = 1, ny
            loopz: do k = 1, nz
                x = x0 + i*dxyz
                y = y0 + j*dxyz
                z = z0 + k*dxyz
                !write (*,*) 'Val = ', fun(x, y, z)
                data_(i, j, k) = fun(x, y, z)
            end do loopz
        end do loopy
    end do loopx

    !Calculate the center values in all 3 directions
    call apply(data_, nx, ny, nz, 0, temp1) !First in x-direction
    call apply(temp1, nx, ny, nz, 1, temp2) !Then in y
    call apply(temp2, nx, ny, nz, 2, a_center) !Then in z

    !call apply(data_, nx, ny, nz, 0, a_center)

    !Now, we test our reconstructed data against the original data
    loopx5: do i = 1, nx
        loopy5: do j = 1, ny
            loopz5: do k = 1, nz
                write (*,*) a_center(i, j, k), data_(i, j, k)
            end do loopz5
        end do loopy5
    end do loopx5



contains

!Dummy function
real*8 function fun(x, y, z)
    implicit none
    real*8, intent(in) :: x, y, z

    fun = 2.*x + 3. * y + 4. * z
end function fun

!Apply function as in weno.c++
subroutine apply(data, nx, ny, nz, dirn, a_center_xyz)
    implicit none
    
    !Declare variables
    integer, intent(in) :: nx, ny, nz, dirn
    real*8, dimension(nx, ny, nz), intent(in) :: data
    real*8, dimension(nx, ny, nz), intent(out) :: a_center_xyz 

    real*8, dimension(nx, ny, nz) :: temp_xyz 

    real*8, dimension(5) :: A = 0
    

    real*8 :: A0, A1, A2, A3, A4, beta1, beta2, beta3, wbarplus1, wbarplus2, wbarplus3, iwbarplussum, wplus1, wplus2, wplus3

    integer :: i, j, k, p, q
    integer, dimension(5,3) :: ijk

    real*8, dimension(3,6) :: beta_shu
    real*8, dimension(3,5) ::  weno_coeffs_de_avg

    
    !Define beta_shu (same for all in WHAM) and WENO coeffs (for de-averaging here)
    beta_shu = reshape((/4.0/3.0, 4.0/3.0, 10.0/3.0,&
                        -19.0/3.0, -13.0/3.0, -31.0/3.0,&
                        25.0/3.0, 13.0/3.0, 25.0/3.0,&
                        11.0/3.0, 5.0/3.0, 11.0/3.0,&
                        -31.0/3.0, -13.0/3.0, -19.0/3.0,&
                        10.0/3.0, 4.0/3.0, 4.0/3.0/), shape(beta_shu))
    
    weno_coeffs_de_avg = reshape((/ -1.0/24.0, 0.,0.,&
                                    1.0/12.0, -1.0/24.0, 0.,&
                                    23.0/24.0, 13.0/12.0, 23.0/24.0,&
                                    0.,-1.0/24.0,1.0/12.0,&
                                    0.,0.,-1.0/24.0&
    /), shape(weno_coeffs_de_avg))

    temp_xyz = 0
    
    loopz: do k = GRHydro_stencil, nz - GRHydro_stencil+1 
        loopy: do j = GRHydro_stencil, ny - GRHydro_stencil+1 
            loopx: do i = GRHydro_stencil, nx - GRHydro_stencil+1 
                !This is in place of the definition of ijk[5] in weno.C++. F90 does not have ternary operators, so I had to do this in a complicated way
                !If there's a better way to do this, please let me know
                select case (dirn)
                case (0) !Solve by x-direction
                    ijk = reshape((/i-2, i-1, i, i+1, i+2, j, j, j, j, j, k, k, k, k, k/), shape(ijk))
                case (1) !Solve by y-direction
                    ijk = reshape((/j, j, j, j, j, i-2, i-1, i, i+1, i+2, k, k, k, k, k/), shape(ijk))
                case (2) !Solve by z-direction
                    ijk = reshape((/k, k, k, k, k, j, j, j, j, j, i-2, i-1, i, i+1, i+2/), shape(ijk))
                end select

                !Again, Cannot define functions inside a function in F90, so yet another convoluted way to do something simple
                A0 = data(ijk(1, 1), ijk(1,2), ijk(1, 3))
                A1 = data(ijk(2, 1), ijk(2,2), ijk(2, 3))
                A2 = data(ijk(3, 1), ijk(3,2), ijk(3, 3))
                A3 = data(ijk(4, 1), ijk(4,2), ijk(4, 3))
                A4 = data(ijk(5, 1), ijk(5,2), ijk(5, 3))

                A = reshape((/A0, A1, A2, A3, A4/), shape(A))

                

                !Calculate smoothness factors
                beta1  = beta_shu(1,1)*SQR(A(1))&
                  + beta_shu(1,2)*A(1)*A(2)&
                  + beta_shu(1,3)*SQR(A(2))&
                  + beta_shu(1,4)*A(1)*A(3)&
                  + beta_shu(1,5)*A(2)*A(3)&
                  + beta_shu(1,6)*SQR(A(3))
                
                beta2  = beta_shu(2,1)*SQR(A(2))&
                  + beta_shu(2,2)*A(2)*A(3)&
                  + beta_shu(2,3)*SQR(A(3))&
                  + beta_shu(2,4)*A(2)*A(4)&
                  + beta_shu(2,5)*A(3)*A(4)&
                  + beta_shu(2,6)*SQR(A(4))
                beta3  = beta_shu(3,1)*SQR(A(3))&
                  + beta_shu(3,2)*A(3)*A(4)&
                  + beta_shu(3,3)*SQR(A(4))&
                  + beta_shu(3,4)*A(3)*A(5)&
                  + beta_shu(3,5)*A(4)*A(5)&
                  + beta_shu(3,6)*SQR(A(5))
                
                wbarplus1 = -9.0/80.0 / SQR(weno_eps + beta1)
                wbarplus2 = 49.0/40.0 / SQR(weno_eps + beta2)
                wbarplus3 = -9.0/80.0 / SQR(weno_eps + beta3)
           
                iwbarplussum = 1.0 / (wbarplus1 + wbarplus2 + wbarplus3)
           
                wplus1 = wbarplus1 * iwbarplussum
                wplus2 = wbarplus2 * iwbarplussum
                wplus3 = wbarplus3 * iwbarplussum
                

                !Calculate the reconstruction
                temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) = 0

                looprecon: do p = 1, 5
                    temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) = &
                    temp_xyz(ijk(3, 1), ijk(3,2), ijk(3, 3)) + ((wplus1 * weno_coeffs_de_avg(1,p) + &
                    wplus2 * weno_coeffs_de_avg(2,p) + wplus3 * weno_coeffs_de_avg(3,p)) * A(p))
                end do looprecon
            end do loopx
        end do loopy
    end do loopz
    a_center_xyz = temp_xyz
    
end subroutine apply

!Gives square of real*8 number
real*8 function SQR(x)
    implicit none
    real*8, intent(in)::x
    SQR = x*x
end function SQR


end program weno