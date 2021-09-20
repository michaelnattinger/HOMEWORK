! ************************************************************************
! Filename : neogrowth_parallel.f90
!
! Author : Philip Coyle
!
! Date Created : September 7th, 2021
!
! Description : This program will use dynamic programming techniques to solve
! a simple neoclassical growth model with a two state markov productivity shock.
!
! Routine:
! cd /Users/philipcoyle/Documents/School/University_of_Wisconsin/ThirdYear/Fall_2021/TA\ -\ Computation/ProblemSets/PS1/Fortran
! gfortran -fopenmp -O2 -o neogrowth_parallel neogrowth_parallel.f90
! ifort -qopenmp -O2 -o neogrowth_parallel neogrowth_parallel.f90
! ./neogrowth_parallel
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
double precision, parameter :: 				cBET 				= 0.99d0
double precision, parameter :: 				cTHETA 			= 0.36d0
double precision, parameter :: 				cDEL 				= 0.025d0

! Model Probabilities
double precision, parameter :: 				Pgg 				= 0.977d0
double precision, parameter :: 				Pbg 				= 1d0 - Pgg
double precision, parameter :: 				Pbb 				= 0.926d0
double precision, parameter :: 				Pgb 				= 1 - Pbb
double precision	 				  :: 				Pr(2)


! Tolerance level for convergence and max itations
double precision, parameter :: 				tol			 		= 1d-4!
integer, parameter 					:: 				max_it 			= 10000
integer 										::				it		 			= 1 !itation counter
integer 										:: 				converged		= 0


! -----------------------------------------------------------------------
! ****************************GRID SET**********************************
! -----------------------------------------------------------------------
! Set up for discritizing the state space (Capital Grid)
integer						 				  :: 				i_k, i_kpr
integer, parameter 				  :: 				n_k 				= 1000
double precision 						:: 				grid_k(n_k)
double precision, parameter :: 				min_k 			= 0.01d0 !1d0
double precision, parameter :: 				max_k 			= 45d0 !75d0
double precision, parameter :: 				step_k 			= (max_k - min_k)/(dble(n_k) - 1d0)
double precision					  :: 				k_today
double precision					  :: 				k_tomorrow

! Set up for discritizing the state space (Productivity Grid)
integer						 				  :: 				i_Z
double precision, parameter :: 				Zg 					= 1.25d0
double precision, parameter :: 				Zb 					= 0.2d0
integer, parameter 				  :: 				n_Z 				= 2
double precision 						:: 				grid_Z(n_Z)
double precision					  :: 				Z_today

integer 										:: 				i_state
integer, parameter  				:: 				n_state 		= n_k*n_Z ! Writing grid in one whole loop

! Global variables for Dynamic Progamming
double precision 						:: 				c_today
double precision 						:: 				c_today_temp
double precision 						:: 				y_today
double precision 						:: 				k_tomorrow_max
double precision 						:: 				v_today
double precision 						:: 				v_today_temp
double precision 						:: 				v_tomorrow


! Allocating space for Policy Functions
double precision 						:: 				pf_c(n_k, n_Z)
double precision 						:: 				pf_k(n_k, n_Z)
double precision 						:: 				pf_v(n_k, n_Z)

integer 										::			  i_stat
integer 										::			  iMaxThreads

end module params_grid


! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : neogrowth
!
! Description : This program will use dynamic programming techniques to solve
! a simple neoclassical growth model with a two state markov productivity shock
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

program neogrowth

use params_grid
use omp_lib

implicit none

! Begin Computational Timer
integer 										::				beginning, end, rate

call system_clock(beginning, rate)

! Initialize Grids and Policy Functions
call housekeeping()

! Do Value Function Iteration
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

end program neogrowth

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
use omp_lib

implicit none

! Discretizing the state space (capital)
do i_k = 1,n_k
	grid_k(i_k) = min_k + (dble(i_k) - 1d0)*step_k
end do

! Discretizing the state space (productivity)
do i_Z = 1,n_Z
	if (i_Z == 1) then
		grid_Z(i_Z) = Zg
	else
		grid_Z(i_Z) = Zb
	end if
end do

! Setting up Policy Function guesses
do i_k = 1,n_k
	do i_Z = 1,n_Z
		pf_c(i_k, i_Z) 			    = 0d0
		pf_k(i_k, i_Z) 		    	= 0d0
		pf_v(i_k, i_Z)    	    = 0d0
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
use omp_lib

implicit none
! allocating space for policy function updates
double precision 						:: 				pf_c_up(n_k, n_Z)
double precision 						:: 				pf_k_up(n_k, n_Z)
double precision 						:: 				pf_v_up(n_k, n_Z)

double precision 						:: 				diff_c
double precision 						:: 				diff_k
double precision 						:: 				diff_v
double precision 						:: 				max_diff


converged = 0
it = 1

! Begin Dynamic Programming Algo
do while (converged == 0 .and. it < max_it)

	!iMaxThreads = omp_get_max_threads()
	!call omp_set_num_threads(iMaxThreads) ! iMaxThreads = 2*(num of cores)


	! When parallelizing over one grid, need to include state vars in private.
	! Break OMP line up for gfortran compiler
	!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i_state, i_Z, i_K, Z_today, Pr, k_today, k_tomorrow, y_today, c_today_temp, &
	!$OMP 																 v_tomorrow, v_today_temp, v_today, k_tomorrow_max, c_today)
	!$OMP DO
	do i_state = 1,n_state

		! ***********************************
		! Define the today's state variables
		! ***********************************
		i_Z         = MOD(i_state,n_Z) + 1
		i_k         = MOD(CEILING(dble(i_state)/dble(n_Z)),n_k) + 1

		Z_today = grid_Z(i_Z)
		k_today = grid_k(i_k)

		! *********************
		! Define today's shock
		! *********************
		if (i_Z == 1) then
			Pr(1) = Pgg
			Pr(2) = Pbg
		else
			Pr(1) = Pgb
			Pr(2) = Pbb
		end if

		! ******************************************************
		! Solve for the optimal consumption / capital investment
		! ******************************************************
		v_today = -1d10
		y_today = Z_today*k_today**(cTHETA)
		do i_kpr = 1,n_k
			k_tomorrow = grid_k(i_kpr)

			c_today_temp = y_today + (1-cDEL)*k_today - k_tomorrow
			v_tomorrow = Pr(1)*pf_v(i_kpr,1) + Pr(2)*pf_v(i_kpr,2)


			c_today_temp = max(0d0,c_today_temp)
			v_today_temp = log(c_today_temp) + cBET*v_tomorrow

			if (v_today_temp > v_today) then
				v_today = v_today_temp
				k_tomorrow_max = k_tomorrow
			end if

		end do

		k_tomorrow = k_tomorrow_max
		c_today = y_today + (1-cDEL)*k_today - k_tomorrow

   	! *******************************
	  ! ****Update Policy Functions****
	  ! *******************************
	  pf_c_up(i_k, i_Z) = c_today
	  pf_k_up(i_k, i_Z) = k_tomorrow
	  pf_v_up(i_k, i_Z) = v_today
	end do
	!$OMP END DO
	!$OMP END PARALLEL


	! Find the difference between the policy functions and updates
	diff_c  = maxval(abs(pf_c - pf_c_up))
	diff_k  = maxval(abs(pf_k - pf_k_up))
	diff_v  = maxval(abs(pf_v - pf_v_up))

	max_diff = diff_c + diff_k + diff_v

	if (mod(it,50) == 0) then
		write(*,*) ""
		write(*,*) "********************************************"
		write(*,*) "At itation = ", it
		write(*,*) "Max Difference = ", max_diff
		write(*,*) "********************************************"
	 end if

	if (max_diff < tol) then
		converged = 1
		write(*,*) ""
		write(*,*) "********************************************"
		write(*,*) "At itation = ", it
		write(*,*) "Max Difference = ", max_diff
		write(*,*) "********************************************"
	end if

	it = it+1

	pf_c 		= pf_c_up
	pf_k 		= pf_k_up
	pf_v		= pf_v_up
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
use omp_lib

implicit none

write(*,*) ""
write (*,*) "Writing PFs to DAT file"
open(unit = 2, file = 'pfs_neogrowth.dat', status = 'replace', action = 'write', iostat = i_stat)
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x)

do i_Z = 1, n_Z
	do i_k = 1,n_k
	     write(2,200) grid_Z(i_Z),grid_k(i_k), pf_c(i_k, i_Z), pf_k(i_k, i_Z), pf_v(i_k, i_Z)
	end do
end do

return

end subroutine coda
