program driver
	
	use forward_diff, n => nx, fp => forward_problem !! Alias to a local variable
	use forward_tgt, n_useless => nx, fp_useless => forward_problem

	implicit none 

	real(8), dimension(n+1) :: xx=0., xx_tlm =0., xxb =0., xx_fd, accuracy
	real(8) :: V=0., Vb = 1., V_forward, Vd = 0.
	real(8), parameter :: eps = 1.d-7
	integer :: ii

	call fp(xx,V)

	!! Forward run
	V_forward = V


	!! Adjoint run
	xx = 0.
	V = 0.
	call forward_problem_b(xx,xxb,V,Vb)

	print *, "#    Reverse    FD    Tangent    Relative accuracy"
	!! Finite differences and Tangent Linear Model
	do ii = 1, n+1

		xx = 0.
		V = 0.
		xx_tlm = 0.
		xx_tlm(ii) = 1.

		!! TLM
		call forward_problem_d(xx,xx_tlm,V,Vd)


		!! FD
		xx = 0.
		V = 0.
		xx(ii) = eps
		call fp(xx,V)
		xx_fd(ii) =  (V - V_forward)/eps

        if ( xx_fd(ii).NE. 0. ) then
            accuracy(ii) = 1.d0 - Vd/xx_fd(ii)
        else
            accuracy(ii) = 0.
        end if
        print *, ii, "    ", xxb(ii), "    ", xx_fd(ii),"    ", Vd,"    ", accuracy(ii)
	end do
end program driver

