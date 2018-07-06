! This source file (a1.f90) was written by Atsushi MOGI (s-182630) for problem 2 A-1,
! and was downloaded from https://github.com/mogi-atsushi/matrix_report/blob/master/a1.f90
! In order to compile this source file (a1.f90), run as follows:
! gfortran a1.f90 -L ~/lapack-3.8.0/ -llapack -lrefblas

program matrix_report_a1
  implicit none
  integer :: nx = 20, ny = 20, i, j, k, n, info
  integer, allocatable :: f(:,:,:), ipiv(:)
  real(8), allocatable :: a(:,:), b(:), u(:), ux(:,:), uy(:,:)
  real(8) :: lambda = 1.0_8, mu = 1.0_8, dx, dy, c = -1.0_8

  dx = 1.0_8/nx
  dy = 1.0_8/ny

  allocate(f(nx+1, ny+1, 2))

  ! Set f(i,j,k)
  n = 0
  do i = 1, nx + 1
     do j = 1, ny + 1
        do k = 1, 2
           n = n + 1
           f(i,j,k) = n
        end do
     end do
  end do

  allocate(a(n, n))
  allocate(b(n))

  do i = 1, n
     b(i) = 0.0_8
     do j = 1, n
        a(i,j) = 0.0_8
     end do
  end do

  ! Equations of motion
  do i = 2, nx
     do j = 2, ny
        a(f(i,j,1), f(i-1,j,1)) = (lambda + 2.0_8*mu)/(dx**2.0_8)
        a(f(i,j,1), f(i+1,j,1)) = (lambda + 2.0_8*mu)/(dx**2.0_8)
        a(f(i,j,1), f(i,j,1)) = -(2.0_8*lambda + 4.0_8*mu)/(dx**2) -(2.0_8*mu)/(dy**2.0_8)
        a(f(i,j,1), f(i,j-1,1)) = mu/(dy**2.0_8)
        a(f(i,j,1), f(i,j+1,1)) = mu/(dy**2.0_8)
        a(f(i,j,1), f(i-1,j-1,2)) = (lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,1), f(i-1,j+1,2)) = -(lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,1), f(i+1,j-1,2)) = -(lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,1), f(i+1,j+1,2)) = (lambda + mu)/(4.0_8*dx*dy)

        a(f(i,j,2), f(i,j-1,2)) = (lambda + 2.0_8*mu)/(dy**2.0_8)
        a(f(i,j,2), f(i,j+1,2)) = (lambda + 2.0_8*mu)/(dy**2.0_8)
        a(f(i,j,2), f(i,j,2)) = -(2.0_8*lambda + 4.0_8*mu)/(dy**2) -(2*mu)/(dx**2.0_8)
        a(f(i,j,2), f(i-1,j,2)) = mu/(dx**2)
        a(f(i,j,2), f(i+1,j,2)) = mu/(dx**2)
        a(f(i,j,2), f(i-1,j-1,1)) = (lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,2), f(i-1,j+1,1)) = -(lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,2), f(i+1,j-1,1)) = -(lambda + mu)/(4.0_8*dx*dy)
        a(f(i,j,2), f(i+1,j+1,1)) = (lambda + mu)/(4.0_8*dx*dy)
     end do
  end do

  ! Boundary conditions
  do i = 2, nx
     do k = 1, 2
        a(f(i,ny+1,k), f(i,ny+1,k)) = 1.0_8
        a(f(i,ny+1,k), f(i,2,k)) = -1.0_8
        a(f(i,1,k), f(i,1,k)) = 1.0_8
        a(f(i,1,k), f(i,ny,k)) = -1.0_8
     end do
  end do
  do j = 1, ny + 1
     do k = 1, 2
        a(f(1,j,k), f(1,j,k)) = 1.0_8
        a(f(nx+1,j,k), f(nx+1,j,k)) = 1.0_8
     end do
  end do
  do j = 1, ny + 1
     b(f(nx+1,j,1)) = c
  end do

  allocate(ipiv(n))
  allocate(u(n))

  u = b

  deallocate(b)

  ! Solve a u = b
  call dgesv(n, 1, a, n, ipiv, u, n, info)

  deallocate(a)
  deallocate(ipiv)

  u = u * 0.03_8

  ! Output
  do i = 1, nx + 1
     do j = 1, ny + 1
        write(*,*) (i-1)*dx, (j-1)*dy, u(f(i,j,1)), u(f(i,j,2))
     end do
  end do

  deallocate(f)
  deallocate(u)
  stop
end program matrix_report_a1
