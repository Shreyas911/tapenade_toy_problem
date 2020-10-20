program driver
	
	use forward_diff
	implicit none 
	real(8), dimension(nx+1) :: xx=0., xxb=0.
	real(8) :: V=0., Vb = 1.
	call forward_problem(xx,V)
	call forward_problem_b(xx,xxb,V,Vb)
	print *,xxb

end program driver

