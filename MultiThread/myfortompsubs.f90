subroutine finit(N, scl, vec) bind(C)

	use iso_c_binding
	implicit none

	! input variables
	integer(c_int), value, intent(in) :: N
	real(c_double), value, intent(in) :: scl
	real(c_double), dimension(N), intent(inout) :: vec

	! thread variables
	integer(c_int) :: omp_get_num_threads, omp_get_thread_num
	integer(c_int) :: nthds, tid

	! compute variables
	integer(c_int) :: stride, start, stp, ii

!$omp parallel private(nthds, tid, stride, start, stp, ii) shared(N, scl, vec)

		! compute thread variables
		nthds = omp_get_num_threads()
		tid = omp_get_thread_num()

		! compute stride
		stride = ceiling(dble(N)/dble(nthds))

		! compute start and stop
		start = tid*stride + 1
		stp = min((tid+1)*stride,N)

		! print info
		!print*, "id, stride, start, stop = ",tid,stride,start,stp

		! initialize vec
		do ii=start,stp
			vec(ii) = scl
		end do

!$omp end parallel

end subroutine

subroutine fclenshaw(N, coeff, vecin, alpha, temp2, beta, vecout, scl, temp1) bind(C)

	use iso_c_binding
	implicit none

	! input variables
	integer(c_int), value, intent(in) :: N
	real(c_double), value, intent(in) :: coeff, alpha, beta, scl
	real(c_double), dimension(N), intent(in) :: vecin, temp2
	real(c_double), dimension(N), intent(inout) :: vecout, temp1

	! thread variables
	integer(c_int) :: omp_get_num_threads, omp_get_thread_num
	integer(c_int) :: nthds, tid

	! compute variables
	integer(c_int) :: stride, start, stp, ii
	real(c_double) :: temp

!$omp parallel private(nthds, tid, stride, start, stp, ii, temp) shared(N, coeff, vecin, alpha, temp2, beta, vecout, scl, temp1)

		! compute thread variables
		nthds = omp_get_num_threads()
		tid = omp_get_thread_num()

		! compute stride
		stride = ceiling(dble(N)/dble(nthds))

		! compute start and stop
		start = tid*stride + 1
		stp = min((tid+1)*stride,N)

		! print info
		!print*, "id, stride, start, stop = ",tid,stride,start,stp

		! clenshaw vec
		do ii=start,stp
			temp = vecout(ii)
			vecout(ii) = coeff*vecin(ii) + alpha*temp2(ii) + beta*vecout(ii) + scl*temp1(ii)
			temp1(ii) = temp
		end do

!$omp end parallel
end subroutine

subroutine fmult(Dims, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax, specrad, vecin, vecout) bind(C)

	use iso_c_binding
	implicit none

	! input variables
	integer(c_int), value, intent(in) :: Dims, Nx, Ny, Nz
	real(c_double), value, intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax, specrad
	real(c_double), dimension(Nx*NY*Nz), intent(in) :: vecin
	real(c_double), dimension(Nx*Ny*Nz), intent(inout) :: vecout

	! thread variables
	integer(c_int) :: omp_get_num_threads, omp_get_thread_num
	integer(c_int) :: nthds, tid

	! compute variables
	integer(c_int) :: N, stride, start, stp, ii
	integer(c_int) :: xpos, ypos, zpos
	real(c_double) :: hx, hy, hz

!$omp parallel private(nthds, tid, N, stride, start, stp, ii, xpos, ypos, zpos, hx, hy, hz), & 
!$omp& shared(Dims, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax, specrad, vecin, vecout)

		! compute thread variables
		nthds = omp_get_num_threads()
		tid = omp_get_thread_num()

		! compute problem size 
		N = Nx*Ny*Nz

		! compute stride
		stride = ceiling(dble(N)/dble(nthds))

		! compute start and stop
		start = tid*stride + 1
		stp = min((tid+1)*stride,N)

		! print info
		!print*, "id, stride, start, stop = ",tid,stride,start,stp

		! compute directional differences
		hx = dble(xmax - xmin)/dble(Nx + 1)
		hy = dble(ymax - ymin)/dble(Ny + 1)
		hz = dble(zmax - zmin)/dble(Nz + 1)

		! loop for multiplication
		do ii=start,stp
			! compute xpos, ypos, zpos
			xpos = mod(ii,Nx)
			if(xpos == 0)then 
				xpos = Nx 
			end if
			ypos = mod((ii-xpos)/Nx,Ny)+1
			zpos = ((ii-xpos)/Nx-ypos+1)/Ny+1
			!print*, "id, xpos, ypos, zpos = ",tid,xpos,ypos,zpos

			! negative x laplacian
			if(xpos == 1 .AND. ii<=N)then
				vecout(ii) = (2d0*vecin(ii) - vecin(ii+1))/hx/hx/specrad
			else if(xpos == Nx .AND. ii<=N)then
				vecout(ii) = (2d0*vecin(ii) - vecin(ii-1))/hx/hx/specrad
			else if(ii <= N)then
				vecout(ii) = (2d0*vecin(ii) - vecin(ii+1) - vecin(ii-1))/hx/hx/specrad
			end if

			! negative y laplacian
			if(Dims == 2)then
				if(ypos == 1 .AND. ii<=N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii+Nx))/hy/hy/specrad
				else if(ypos == Ny .AND. ii<=N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii-Nx))/hy/hy/specrad
				else if(ii <= N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii+Nx) - vecin(ii-Nx))/hy/hy/specrad
				end if
			end if

			! negative z laplacian
			if(Dims == 3)then
				if(zpos == 1 .AND. ii<=N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii+Nx*Ny))/hz/hz/specrad
				else if(zpos == Nz .AND. ii<=N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii-Nx*Ny))/hz/hz/specrad
				else if(ii <= N)then
					vecout(ii) = vecout(ii) + (2d0*vecin(ii) - vecin(ii+Nx*Ny) - vecin(ii-Nx*Ny))/hz/hz/specrad
				end if
			end if
		end do

!$omp end parallel
end subroutine
