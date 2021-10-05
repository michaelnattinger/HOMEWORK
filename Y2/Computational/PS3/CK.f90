! ************************************************************************
! Filename : neogrowth_deterministic.f90
! Author : Michael Nattinger
!
! Date Created : September 7th, 2021
!
! Description : This program will use dynamic programming techniques to solve
! a simple neoclassical growth model with a two state markov productivity shock.
!
! Notes:
! I run my fortan code from the terminal. I included the commands I need to run
! in the terminal.
! First, I chage into the appropriate directory.
! Second I compile the code. gfortran is my compiler. I will go over what all the dashes mean when we meet.
! Finally, I run the compiled (machine readbale) code. This is running the program.
! You should update the path -- the part that is "cd ..." -- to where you store the file.
!
! Routine:
! cd /Users/philipcoyle/Documents/School/University_of_Wisconsin/ThirdYear/Fall_2021/TA\ -\ Computation/ProblemSets/PS1/Fortran
! gfortran -O2 -o neogrowth_deterministic neogrowth_deterministic.f90
! ./neogrowth_deterministic
! ************************************************************************

! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! module : params_grid
!
! Description : This module will form the foudation for our program. In it
! we will allocate space for all paramaters used in this program and set up the
! grids to create a discritized state space
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

module params_grid

implicit none

! -----------------------------------------------------------------------
! *******************DECLARATION OF PARAMETERS AND VARIABLES*************
! -----------------------------------------------------------------------
! Model Parameters
integer, parameter :: N_lifetime = 66 ! number of periods our agents live for
double precision, parameter :: n_popg = 0.011d0
integer, parameter :: J_retire = 46 ! retirement age
double precision, parameter :: cTHETA = 0.11d0
double precision, parameter :: cGAMMA = 0.42d0
double precision, parameter :: cSIGMA = 2d0
double precision, parameter :: cDELTA = 0.06d0
double precision, parameter :: 				cBET 				= 0.97d0 	! Discount Factor
double precision, parameter :: 				cALPHA 			= 0.36d0	! Capital Share
double precision, parameter :: cERGO_g = 0.2037
double precision, parameter :: cERGO_b = 1-cERGO_g ! ergodic distribution
double precision, parameter :: 				cZh 				= 3.0d0	! High productivity
double precision, parameter :: 				cZl  				= 0.5d0	! Low productivity
double precision, parameter :: 				cPhh 				= 0.9261d0	! Probability of good -> good transition
double precision, parameter :: 				cPhl 				= 1-cPhh	! Probability of good -> bad transition
double precision, parameter :: 				cPll 				= 0.9811d0	! Probability of bad -> good transition
double precision, parameter :: 				cPlh 				= 1-cPll	! Probability of bad -> bad transition
double precision, parameter :: wage = 1.05
double precision, parameter :: rental = 0.05
double precision, parameter :: benefits = 0.2


! Tolerance level for convergence and max itations
double precision, parameter :: 				tol			 		= 1d-6 	! Convergence Tolerance
double precision, parameter :: 				q_tol			 		= 1d-3 	! Convergence Tolerance for a market clearing (q loop tol)
double precision, parameter :: 				pmf_tol			 		= 1d-9 	! Convergence Tolerance for pmf
integer, parameter 					:: 				max_it 			= 10000 ! Maximum Number of Iterations
integer 										::				it		 			= 1 		! Itation Counter
integer 										:: 				converged		= 0			! Dummy for VFI Convergence
integer :: q_converged = 0 ! has q converged?
integer :: pmf_converged = 0
integer, parameter :: n_samples = 5 ! must be >3

! -----------------------------------------------------------------------
! ****************************GRID SET**********************************
! -----------------------------------------------------------------------
! Set up for discritizing the state space 
integer				  :: 				i_a, i_apr									! Ieration Counters for k_today and k_tomorrow grid
integer, parameter 				  :: 				n_a 				= 5000																! Size of k grid
double precision 						:: 				grid_a(n_a)																				! Allocate Space for k grid
double precision, parameter :: 				min_a 			= 0d0															! Minimum of k grid
double precision, parameter :: 				max_a 			= 75d0																! Maximum of k grid
double precision, parameter :: 				step_a 			= (max_a - min_a)/(dble(n_a) - 1d0) 	! Step of k grid
double precision					  :: 				a_today
double precision					  :: 				a_tomorrow
integer :: i_z ! 1 good state, 2 bad state

! Global variables for Dynamic Progamming
double precision 						:: 				c_today
double precision 						:: 				c_today_temp
double precision 						:: 				y_today
double precision 						:: 				a_tomorrow_max
double precision 						:: 				v_today
double precision 						:: 				v_today_temp
double precision 						:: 				v_tomorrow


! Allocating space for Policy Functions
double precision 						:: 				pf_c(n_a,N_lifetime) ! Allocate Space for Conumption Policy Function ! no _b means it is implicitly _g
double precision 						:: 				pf_a(n_a,N_lifetime) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v(n_a,N_lifetime) ! Allocate Space for Value Function
double precision 						:: 				pf_l(n_a,N_lifetime) ! Allocate Space for Value Function

double precision 						:: 				pf_c_b(n_a,N_lifetime) ! Allocate Space for Conumption Policy Function
double precision 						:: 				pf_a_b(n_a,N_lifetime) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v_b(n_a,N_lifetime) ! Allocate Space for Value Function
double precision 						:: 				pf_l_b(n_a,N_lifetime) ! Allocate Space for Value Function

double precision :: pmf_init(n_a,2,N_lifetime)
double precision :: pmf(n_a,2,N_lifetime)
integer						:: 				pf_i_apr(n_a,2,N_lifetime)
!double precision 						:: pmf_init_flat(n_a,2,N_lifetime)
integer 										::			  i_stat   ! Used for writing data after program has run.
integer :: i_age 
! Price
double precision :: q_guess
end module params_grid


! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : huggett
!
! Description : Implements Conesa and Krueger (1993, JEDC) - Michael Nattinger
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

program huggett

use params_grid


implicit none

! Begin Computational Timer
integer 										::				beginning, end, rate

call system_clock(beginning, rate)

! Initialize Grids and Policy Functions
call housekeeping()

! Do Value Function Iteration & Search for Fixed Point
call bellman()

call system_clock(end)
write(*,*) ""
write(*,*) "******************************************************"
write(*,*) "Total elapsed time = ", real(end - beginning) / real(rate)," seconds"
write(*,*) "******************************************************"

! Write results
call coda()

write(*,*) ""
write(*,*) "**************************************"
write(*,*) "************END OF PROGRAM************"
write(*,*) "**************************************"
write(*,*) ""

end program huggett

! ************************************************************************
! ************************************************************************
! ************************************************************************
! **************************** SUBROUTINES *******************************
! ************************************************************************
! ************************************************************************
! ************************************************************************


! ************************************************************************


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : housekeeping
!
! description : Initializes Grids and Policy Functions
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine housekeeping()

use params_grid


implicit none

! Discretizing the state space (capital)
do i_a = 1,n_a
	grid_a(i_a) = min_a + (dble(i_a) - 1d0)*step_a
end do

! Setting up Policy Function guesses
do i_a = 1,n_a
do i_age = 1,N_lifetime
	pf_c(i_a,i_age) 			    = 0d0
	pf_a(i_a,i_age) 		    	= 0d0
	pf_v(i_a,i_age)    	      = 0d0
	pf_l(i_a,i_age)    	      = 0d0

	pf_c_b(i_a,i_age) 			    = 0d0
	pf_a_b(i_a,i_age) 		    	= 0d0
	pf_v_b(i_a,i_age)    	      = 0d0
	pf_l_b(i_a,i_age)    	      = 0d0

	pmf_init(i_a,1,i_age) = 0d0
	pmf_init(i_a,2,i_age) = 0d0
	pf_i_apr(i_a,1,i_age) = 0
	pf_i_apr(i_a,2,i_age) = 0
	!pmf_init_flat(i_a,1) = 1d0/(dble(n_a) * 2d0)
	!pmf_init_flat(i_a,2) = 1d0/(dble(n_a) * 2d0)
	end do
end do

return

end subroutine housekeeping




! ************************************************************************


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : bellman
!
! description : Solves the dynamic programming problem for policy and value
! functions.
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine bellman()

use params_grid


implicit none
! allocating space for policy function updates
double precision 						:: 				pf_c_up(n_a,N_lifetime)
double precision 						:: 				pf_a_up(n_a,N_lifetime)
double precision 						:: 				pf_v_up(n_a,N_lifetime)
double precision 						:: 				pf_l_up(n_a,N_lifetime)
double precision 						:: 				pf_c_up_b(n_a,N_lifetime)
double precision 						:: 				pf_a_up_b(n_a,N_lifetime)
double precision 						:: 				pf_v_up_b(n_a,N_lifetime)
double precision 						:: 				pf_l_up_b(n_a,N_lifetime)

double precision 						:: 				diff_c
double precision 						:: 				diff_a
double precision 						:: 				diff_v
double precision 						:: 				diff_l
double precision 						:: 				diff_c_b
double precision 						:: 				diff_a_b
double precision 						:: 				diff_v_b
double precision 						:: 				diff_l_b
double precision 						:: 				pmf_diff
double precision 						:: 				max_diff
double precision :: l_today_temp
double precision :: l_today_max
double precision :: eff_w
double precision :: Ph
double precision :: Pl
double precision :: sum_pmf
integer :: apr_floor
integer :: next_e
integer :: next_u
integer :: decl
!integer :: n_samples = 10 ! number of samples for iterated grid search
integer :: i_sample
integer :: i_min(1)
integer :: i_mid 
integer :: n_mid = (n_samples+1)/2
integer :: age 
double precision :: apr_sum

double precision :: pmf_prime(n_a,2)
double precision :: q_min = cBET!0.5!cBET!0.99!0.99!cBET - 0.1!0.4d0
double precision :: q_max = 1.0!2.0!0.5!1.0!cBET + 0.05!cBET
double precision :: q_grid(n_samples-2)
double precision :: apr_mid = 1d0
double precision :: age_prod( J_retire-1)

	  age_prod(1) = 0.59923239  ! I do not know how to read data into fortran yet
	  age_prod(2) = 0.63885106 
	  age_prod(3) = 0.67846973 
      age_prod(4) = 0.71808840 
      age_prod(5) = 0.75699959 
      age_prod(6) = 0.79591079 
      age_prod(7) = 0.83482198 
      age_prod(8) = 0.87373318 
      age_prod(9) = 0.91264437 
      age_prod(10) = 0.95155556 
      age_prod(11) = 0.99046676 
      age_prod(12) = 0.99872065 
      age_prod(13) = 1.0069745 
      age_prod(14) =  1.0152284 
      age_prod(15) =  1.0234823 
      age_prod(16) =  1.0317362 
      age_prod(17) =  1.0399901 
      age_prod(18) =  1.0482440 
      age_prod(19) =  1.0564979 
      age_prod(20) =  1.0647518 
      age_prod(21) =  1.0730057 
      age_prod(22) =  1.0787834 
      age_prod(23) =   1.0845611 
      age_prod(24) =  1.0903388 
      age_prod(25) =  1.0961165 
      age_prod(26) =  1.1018943 
      age_prod(27) =  1.1076720 
      age_prod(28) =  1.1134497 
      age_prod(29) =  1.1192274 
      age_prod(30) =  1.1250052 
      age_prod(31) =  1.1307829 
      age_prod(32) =  1.1233544 
      age_prod(33) =  1.1159259 
      age_prod(34) =  1.1084974 
      age_prod(35) =  1.1010689 
      age_prod(36) =  1.0936404 
      age_prod(37) =  1.0862119 
      age_prod(38) =  1.0787834 
      age_prod(39) =  1.0713549 
      age_prod(40) =  1.0639264
      age_prod(41) =  1.0519200
      age_prod(42) =  1.0430000
      age_prod(43) =  1.0363000
      age_prod(44) =  1.0200000
      age_prod(45) =  1.0110000


converged = 0
it = 1
! We use backwards induction here.
do i_age = 1,N_lifetime ! = N+1-i_age period agent - need to work backwards here
	age = N_lifetime + 1 - i_age !
	do i_z = 1,2 	! productivity in current age
	apr_floor = 1													! loop over productivity states
	do i_a = 1,n_a ! Loop over all possible capital states	
		! ***********************************
		! Define the today's state variables
		! ***********************************
		a_today = grid_a(i_a)

		! ******************************************************
		! Solve for the optimal consumption / capital investment
		! ******************************************************
		v_today = -1d10 ! Set a very low bound for value
		! income if working age
		if (age >= J_retire) then 
			y_today = a_today * (1 + rental) + benefits
		end if
		if (age == N_lifetime) then ! we eat our budget, eskeddit
			c_today = y_today
			v_today = ((c_today) ** ((1-cSIGMA)*cGAMMA))/(1-cSIGMA)
			a_tomorrow = 0
		else if (age>=J_retire) then
			if (i_z < 1.5) then ! then I need to calculate opt, otherwise opt is same as i_z = 1 case
			i_apr = apr_floor
			decl = 0
			do while (decl == 0 .and. i_apr<=n_a)
				a_tomorrow = grid_a(i_apr)
				c_today_temp = y_today - a_tomorrow
				if (c_today_temp>0d0) then 
					v_tomorrow = cPhh*pf_v(i_apr,age+1) + cPhl*pf_v_b(i_apr,age+1) 
					v_today_temp = ((c_today_temp)**(cGAMMA)*(1-l_today_temp)**(1-cGAMMA) )**(1-cSIGMA)/(1-cSIGMA) + cBET*v_tomorrow
					if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
						v_today = v_today_temp
						a_tomorrow_max = a_tomorrow
						apr_floor = i_apr
						pf_i_apr(i_a,i_z,age) = apr_floor
					else
						decl = 1
					end if
				end if
				i_apr = i_apr+1
			end do
			a_tomorrow = a_tomorrow_max
			c_today = y_today - a_tomorrow
			l_today_max = 0d0
			else ! 
			a_tomorrow = pf_a_up(i_a,age)
			c_today = pf_c_up(i_a,age)
			v_today = pf_v_up(i_a,age)
			pf_i_apr(i_a,i_z,age) = pf_i_apr(i_a,1,age)
			l_today_max = 0d0
			end if
		else	!working years					
		!do i_apr = apr_floor,n_a ! Loop over all possible capital chocies
		i_apr = apr_floor	! keep this at 1 until I verify increasing a' choice
		decl = 0 ! this can stay	
		a_tomorrow_max =0d0
		if (i_z < 1.5) then 
			eff_w = (1-cTHETA)*wage*cZh*age_prod(age)
		else
			eff_w = (1-cTHETA)*wage*cZl*age_prod(age)
		end if
		do while (decl == 0 .and. i_apr<=n_a)				
			a_tomorrow = grid_a(i_apr)
			l_today_temp = (cGAMMA * eff_w - (1-cGAMMA)*((1+rental)*a_today - a_tomorrow))/eff_w!
			if (l_today_temp<0d0) then
			l_today_temp = 0d0
			else if (l_today_temp>1d0) then
			l_today_temp = 1d0
			end if
			c_today_temp = eff_w*l_today_temp + (1+rental)*a_today - a_tomorrow
			if (c_today_temp>0d0) then 	! if this thing is negative then don't check it
			if (i_z <1.5) then
				v_tomorrow = cPhh*pf_v(i_apr,age+1) + cPhl*pf_v_b(i_apr,age+1)  							
			else
				v_tomorrow = cPlh*pf_v(i_apr,age+1) + cPll*pf_v_b(i_apr,age+1)  
			end if
			v_today_temp = ((c_today_temp)**(cGAMMA)*(1-l_today_temp)**(1-cGAMMA) )**(1-cSIGMA)/(1-cSIGMA) + cBET*v_tomorrow

			if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
				v_today = v_today_temp
				a_tomorrow_max = a_tomorrow
				l_today_max = l_today_temp
				apr_floor = i_apr
				pf_i_apr(i_a,i_z,age) = i_apr!apr_floor
			else
				decl = 1
			end if
			end if
			i_apr = i_apr+1
		end do
		a_tomorrow = a_tomorrow_max
		c_today = eff_w*l_today_max + (1+rental)*a_today - a_tomorrow
		end if
   		! *******************************
	  	! ****Update Policy Functions****
	  	! *******************************
	  	if (i_z < 1.5) then 
	  		pf_c_up(i_a,age) = c_today 												
	  		pf_a_up(i_a,age) = a_tomorrow												
	  		pf_v_up(i_a,age) = v_today	
	  		pf_l_up(i_a,age) = l_today_max
	  	else
	  		pf_c_up_b(i_a,age) = c_today 												 
	  		pf_a_up_b(i_a,age) = a_tomorrow												 
	  		pf_v_up_b(i_a,age) = v_today	
	  		pf_l_up_b(i_a,age) = l_today_max
	  	end if											
	end do
	end do
	pf_c 		= pf_c_up
	pf_a 		= pf_a_up
	pf_v			= pf_v_up
	pf_l			= pf_l_up
	pf_c_b 		= pf_c_up_b
	pf_a_b 		= pf_a_up_b
	pf_v_b		= pf_v_up_b
	pf_l_b			= pf_l_up
end do
pmf = pmf_init
pmf(1,1,1) = cERGO_g
pmf(1,2,1) = cERGO_b
! iterate forwards
do age = 1,(N_lifetime - 1)
	do i_a = 1,n_a
		do i_z = 1,2
			if (i_z<1.5) then
				Ph = cPhh
				Pl = cPhl
			else
				Ph = cPlh 
				Pl = cPll
			end if
			if (pmf(i_a,i_z,age)>0d0) then ! iterate forwards; account for pop growth
				pmf(pf_i_apr(i_a,i_z,age),1,age+1) = pmf(pf_i_apr(i_a,i_z,age),1,age+1) + Ph*pmf(i_a,i_z,age)/(1+n_popg)
				pmf(pf_i_apr(i_a,i_z,age),2,age+1) = pmf(pf_i_apr(i_a,i_z,age),2,age+1) + Pl*pmf(i_a,i_z,age)/(1+n_popg)
			end if
		end do
	end do
end do

sum_pmf = 0 ! normalize the pdf by 
do age = 1,(N_lifetime)
	do i_a = 1,n_a
		do i_z = 1,2
		sum_pmf = sum_pmf + pmf(i_a,i_z,age)
		end do
	end do
end do
do age = 1,(N_lifetime)
	do i_a = 1,n_a
		do i_z = 1,2
			pmf(i_a,i_z,age) = pmf(i_a,i_z,age)/sum_pmf
		end do
	end do
end do

return

end subroutine bellman

! ************************************************************************


! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! subroutine : coda
!
! description : Writes results to .dat file
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

subroutine coda()

use params_grid


implicit none

write(*,*) ""
write (*,*) "Writing PFs to DAT file"
open(unit = 2, file = 'pfs_.dat', status = 'replace', action = 'write', iostat = i_stat)
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x) !,f25.15,2x,f25.15,2x

do i_a = 1,n_a
     write(2,200) pf_v(i_a,50), pf_a(i_a,20), pf_a_b(i_a,20), pf_a(i_a,1), pf_a_b(i_a,2), &
      pmf(i_a,1,N_lifetime), pmf(i_a,2,N_lifetime)
     !, pf_i_apr(i_a,1), pf_i_apr(i_a,2)!grid_a(i_a), pf_c(i_a), pf_a(i_a), pf_v(i_a), pmf(i_a,1), pmf(i_a,2), pmf_init_flat(i_a,1), pmf_init_flat(i_a,2)!, pf_i_apr(i_a,1), pf_i_apr(i_a,2)
end do

return

end subroutine coda
