program driver
	
	use forward_diff
	implicit none 
	real(8), dimension(nx+1) :: xx=0., xx_tlm =0., xxb =0., xx_fd, accuracy
	real(8) :: V=0., Vb = 1., V_forward, Vd = 0.
	real(8), parameter :: eps = 1.d-7
	integer :: ii

	call forward_problem(xx,V)

	!! Forward run
	V_forward = V


	!! Adjoint run
	xx = 0.
	V = 0.
	call forward_problem_b(xx,xxb,V,Vb)
	print *,xxb
	
	print *, "#    Reverse    FD    Tangent    Relative accuracy"
	!! Finite differences and Tangent Linear Model
	do ii = 1, nx+1

		xx = 0.
		V = 0.
		xx_tlm = 0.
		xx_tlm(ii) = 1.

		!! TLM
		! call forward_problem_d(xx,xx_tlm,V,Vd)


		!! FD
		xx = 0.
		V = 0.
		xx(ii) = eps
		call forward_problem(xx,V)
		xx_fd(ii) =  (V - V_forward)/eps

        if ( xx_fd(ii).NE. 0. ) then
            accuracy(ii) = 1.d0 - Vd/xx_fd(ii)
        else
            accuracy(ii) = 0.
        end if
        print *, ii, "    ", xxb(ii), "    ", xx_fd(ii),"    ", Vd,"    ", accuracy(ii)
	end do
end program driver

