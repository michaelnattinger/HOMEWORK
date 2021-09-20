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
double precision, parameter :: 				cBET 				= 0.9932d0 	! Discount Factor
double precision, parameter :: 				cALPHA 			= 1.5d0	! Capital Share
double precision, parameter :: 				cE 				= 1.0d0	! Depreciation Rate
double precision, parameter :: 				cU 				= 0.5d0	! Depreciation Rate
double precision, parameter :: 				cPee 				= 0.97d0	! Probability of good -> good transition
double precision, parameter :: 				cPeu 				= 0.03d0	! Probability of good -> bad transition
double precision, parameter :: 				cPue 				= 0.5d0	! Probability of bad -> good transition
double precision, parameter :: 				cPuu 				= 0.5d0	! Probability of bad -> bad transition


! Tolerance level for convergence and max itations
double precision, parameter :: 				tol			 		= 1d-4 	! Convergence Tolerance
double precision, parameter :: 				q_tol			 		= 1d-2 	! Convergence Tolerance for q
double precision, parameter :: 				pmf_tol			 		= 1d-8 	! Convergence Tolerance for pmf
integer, parameter 					:: 				max_it 			= 10000 ! Maximum Number of Iterations
integer 										::				it		 			= 1 		! Itation Counter
integer 										:: 				converged		= 0			! Dummy for VFI Convergence
integer :: q_converged = 0 ! has q converged?
integer :: pmf_converged = 0


! -----------------------------------------------------------------------
! ****************************GRID SET**********************************
! -----------------------------------------------------------------------
! Set up for discritizing the state space 
integer				  :: 				i_a, i_apr									! Ieration Counters for k_today and k_tomorrow grid
integer, parameter 				  :: 				n_a 				= 1000																! Size of k grid
double precision 						:: 				grid_a(n_a)																				! Allocate Space for k grid
double precision, parameter :: 				min_a 			= -2d0															! Minimum of k grid
double precision, parameter :: 				max_a 			= 5d0																! Maximum of k grid
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
double precision 						:: 				pf_c(n_a) ! Allocate Space for Conumption Policy Function ! no _b means it is implicitly _g
double precision 						:: 				pf_a(n_a) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v(n_a) ! Allocate Space for Value Function

double precision 						:: 				pf_c_b(n_a) ! Allocate Space for Conumption Policy Function
double precision 						:: 				pf_a_b(n_a) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v_b(n_a) ! Allocate Space for Value Function

double precision :: pmf_init(n_a,2)
double precision :: pmf(n_a,2)
integer						:: 				pf_i_apr(n_a,2)
double precision 						:: pmf_init_flat(n_a,2)
integer 										::			  i_stat   ! Used for writing data after program has run.

! Price
double precision :: q_guess
end module params_grid


! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : huggett
!
! Description : Implements Huggett (1993, JEDC) - Michael Nattinger
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
	pf_c(i_a) 			    = 0d0
	pf_a(i_a) 		    	= 0d0
	pf_v(i_a)    	      = 0d0

	pf_c_b(i_a) 			    = 0d0
	pf_a_b(i_a) 		    	= 0d0
	pf_v_b(i_a)    	      = 0d0

	pmf_init(i_a,1) = 0d0
	pmf_init(i_a,2) = 0d0
	pf_i_apr(i_a,1) = 0
	pf_i_apr(i_a,2) = 0
	pmf_init_flat(i_a,1) = 1d0/(dble(n_a) * 2d0)
	pmf_init_flat(i_a,2) = 1d0/(dble(n_a) * 2d0)
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
double precision 						:: 				pf_c_up(n_a)
double precision 						:: 				pf_a_up(n_a)
double precision 						:: 				pf_v_up(n_a)
double precision 						:: 				pf_c_up_b(n_a)
double precision 						:: 				pf_a_up_b(n_a)
double precision 						:: 				pf_v_up_b(n_a)

double precision 						:: 				diff_c
double precision 						:: 				diff_a
double precision 						:: 				diff_v
double precision 						:: 				diff_c_b
double precision 						:: 				diff_a_b
double precision 						:: 				diff_v_b
double precision 						:: 				pmf_diff
double precision 						:: 				max_diff
integer :: apr_floor
integer :: next_e
integer :: next_u
double precision :: apr_sum

double precision :: pmf_prime(n_a,2)
double precision :: q_min = cBET - 0.1!0.4d0
double precision :: q_max = cBET + 0.05!cBET



q_guess = (q_min + q_max)/2!0.9d0!cBET**(2)
q_converged = 0
do while (q_converged == 0 )
converged = 0
it = 1
! Begin Dynamic Programming Algo.
! Continue to iterate while VFI has not converged (converged == 0) and iteration counter is less than maximium numnber of iterations (it < max_it)
do while (converged == 0 .and. it < max_it)
	do i_z = 1,2 	
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
		if (i_z < 1.5) then
			y_today = a_today + cE ! income today for employed
		else
			y_today = a_today + cU
		end if							
		do i_apr = apr_floor,n_a ! Loop over all possible capital chocies					
			a_tomorrow = grid_a(i_apr)

			! some values are "temp" values because we are searching exaustively for the capital/consumption choice that maximizes value
			c_today_temp = y_today - q_guess*a_tomorrow
			if (c_today_temp>0d0) then 	
			if (i_z <1.5) then
				v_tomorrow = cPee*pf_v(i_apr) + cPeu*pf_v_b(i_apr)  							
			else
				v_tomorrow = cPue*pf_v(i_apr) + cPuu*pf_v_b(i_apr)  
			end if
			c_today_temp = max(0d0,c_today_temp)
			v_today_temp = ((c_today_temp)**(1-cALPHA) - 1)/(1-cALPHA) + cBET*v_tomorrow

			if (v_today_temp > v_today) then ! if "temp" value is best so far, record value and capital choice
				v_today = v_today_temp
				a_tomorrow_max = a_tomorrow
				apr_floor = i_apr
				pf_i_apr(i_a,i_z) = apr_floor
			end if
			end if
		end do

		a_tomorrow = a_tomorrow_max
		c_today = y_today  - q_guess*a_tomorrow_max

   		! *******************************
	  	! ****Update Policy Functions****
	  	! *******************************
	  	if (i_z < 1.5) then 
	  		pf_c_up(i_a) = c_today 												
	  		pf_a_up(i_a) = a_tomorrow												
	  		pf_v_up(i_a) = v_today	
	  	else
	  		pf_c_up_b(i_a) = c_today 												 
	  		pf_a_up_b(i_a) = a_tomorrow												 
	  		pf_v_up_b(i_a) = v_today	
	  	end if											
	end do
	end do
	! Find the difference between the policy functions and updates
	diff_c  = maxval(abs(pf_c - pf_c_up))									
	diff_a  = maxval(abs(pf_a - pf_a_up))										
	diff_v  = maxval(abs(pf_v - pf_v_up))										
	diff_c_b  = maxval(abs(pf_c_b - pf_c_up_b))										
	diff_a_b  = maxval(abs(pf_a_b - pf_a_up_b))										
	diff_v_b  = maxval(abs(pf_v_b - pf_v_up_b))							

	!max_diff = diff_c + diff_a + diff_v + diff_c_b + diff_a_b + diff_v_b 	
	max_diff = max(diff_c,diff_c_b)
	!if (mod(it,200) == 0) then
	!	write(*,*) ""
	!	write(*,*) "********************************************"
	!	write(*,*) "At itation = ", it
	!	write(*,*) "Max Difference = ", max_diff
	!	write(*,*) "********************************************"
	!end if

	if (max_diff < tol) then
		converged = 1
		!write(*,*) ""
		!write(*,*) "********************************************"
		!write(*,*) "At itation = ", it
		!write(*,*) "Max Difference = ", max_diff
		!write(*,*) "********************************************"
	end if

	it = it+1

	pf_c 		= pf_c_up
	pf_a 		= pf_a_up
	pf_v			= pf_v_up
	pf_c_b 		= pf_c_up_b
	pf_a_b 		= pf_a_up_b
	pf_v_b		= pf_v_up_b
end do


! perform mu operator 
pmf = pmf_init_flat
pmf_converged = 0
! mu operator here
do while (pmf_converged == 0)
pmf_prime = pmf_init
do i_a = 1,n_a
! explicitly keep reference out of indices
next_e = pf_i_apr(i_a,1)
next_u = pf_i_apr(i_a,2)
pmf_prime(next_e,1) = pmf_prime(next_e,1) + pmf(i_a,1)*cPee
pmf_prime(next_e,2) = pmf_prime(next_e,2) + pmf(i_a,1)*cPeu
pmf_prime(next_u,1) = pmf_prime(next_u,1) + pmf(i_a,2)*cPue
pmf_prime(next_u,2) = pmf_prime(next_u,2) + pmf(i_a,2)*cPuu
end do
pmf_diff = maxval(abs(pmf_prime - pmf))
if (pmf_diff<pmf_tol) then
pmf_converged = 1
else
pmf = pmf_prime
end if
end do


! Check market clearing
apr_sum = 0
do i_a = 1,n_a
apr_sum = apr_sum + pmf(i_a,1)*pf_a(i_a)
apr_sum = apr_sum + pmf(i_a,2)*pf_a_b(i_a)
end do
!apr_sum = 0d0
! check for market clearing q_tol
if (abs(apr_sum) <q_tol) then
q_converged = 1
write(*,*) ""
write(*,*) "********************************************"
write(*,*) "a sum =  ", apr_sum
write(*,*) "final q  = ", q_guess
write(*,*) "********************************************"
else if (apr_sum>0) then ! bond price too low
q_min = q_guess
q_guess = (q_min + q_max)/2
write(*,*) ""
write(*,*) "********************************************"
write(*,*) "a sum =  ", apr_sum
write(*,*) "new q  = ", q_guess
write(*,*) "********************************************"
else 
q_max = q_guess
q_guess = (q_min + q_max)/2
write(*,*) ""
write(*,*) "********************************************"
write(*,*) "a sum =  ", apr_sum
write(*,*) "new q  = ", q_guess
write(*,*) "********************************************"
end if
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
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)

do i_a = 1,n_a
     write(2,200) grid_a(i_a), pf_c(i_a), pf_a(i_a), pf_v(i_a), pmf(i_a,1), pmf(i_a,2), pmf_init_flat(i_a,1), pmf_init_flat(i_a,2)!, pf_i_apr(i_a,1), pf_i_apr(i_a,2)
end do

return

end subroutine coda
