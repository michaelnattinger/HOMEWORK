! ************************************************************************
! 
! Author : Michael Nattinger
!
! Implements transition between two steady states: from with SS to without.
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
double precision :: cTHETA = 0.11d0 !no longer a parameter
double precision :: cGAMMA = 0.42d0
double precision, parameter :: cSIGMA = 2d0
double precision, parameter :: cDELTA = 0.06d0
double precision, parameter :: cBET 				= 0.97d0 	! Discount Factor
double precision, parameter :: cALPHA 			= 0.36d0	! Capital Share
double precision, parameter :: cERGO_g = 0.2037
double precision, parameter :: cERGO_b = 1-cERGO_g ! ergodic distribution
double precision :: cZh 				= 3.0d0	! High productivity
double precision, parameter :: cZl  				= 0.5d0	! Low productivity
double precision, parameter :: cPhh 				= 0.9261d0	! Probability of good -> good transition
double precision, parameter :: cPhl 				= 1-cPhh	! Probability of good -> bad transition
double precision, parameter :: cPll 				= 0.9811d0	! Probability of bad -> good transition
double precision, parameter :: cPlh 				= 1-cPll	! Probability of bad -> bad transition
! still define these here but not as a parameter - these will be edited in the loop as appropriate
double precision :: wage = 1.05
double precision :: rental = 0.05
double precision :: benefits = 0.2
double precision :: working_mass
integer, parameter :: n_exp = 6
double precision :: K_0ss
double precision :: K_Tss
double precision :: L_0ss
double precision :: L_Tss
integer, parameter :: N_transition = 150 ! number of periods over which transition will take place. Old SS at 0, at 1 change.
double precision :: K_0transition(N_transition)
double precision :: K_1transition(N_transition)
double precision :: L_0transition(N_transition)
double precision :: L_1transition(N_transition)


! Tolerance level for convergence and max itations
double precision :: 				tune_K			 		= 0.995d0!0.99d0
double precision :: 				tune_L			 		= 0.995d0!0.99d0
double precision, parameter :: 				a_tol			 		= 1d-4
double precision, parameter :: 				l_tol			 		= 1d-4  
double precision, parameter :: 				tol			 		= 1d-3 	! Convergence Tolerance on capital transition path
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
integer, parameter 				  :: 				n_a 				= 500																! Size of k grid
integer, parameter :: n_z = 2
double precision 						:: 				grid_a(n_a)																				! Allocate Space for k grid
double precision, parameter :: 				min_a 			= 0d0															! Minimum of k grid
double precision, parameter :: 				max_a 			= 60d0																! Maximum of k grid
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
double precision 						:: 				pf_c(n_a,n_z,N_lifetime) ! Allocate Space for Conumption Policy Function ! no _b means it is implicitly _g
double precision 						:: 				pf_a(n_a,n_z,N_lifetime) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v(n_a,n_z,N_lifetime) ! Allocate Space for Value Function
double precision 						:: 				pf_l(n_a,n_z,N_lifetime) ! Allocate Space for Value Function


! Allocate space for transition path functions
double precision :: tran_pf_v(n_a,n_z,N_lifetime,N_transition+1)
double precision :: tran_pf_a(n_a,n_z,N_lifetime,N_transition)
double precision :: tran_pf_c(n_a,n_z,N_lifetime,N_transition)
double precision :: tran_pf_l(n_a,n_z,N_lifetime,N_transition)
double precision :: tran_pmf(n_a,n_z,N_lifetime,N_transition)
integer :: tran_pf_i_apr(n_a,n_z,N_lifetime,N_transition)


double precision :: pmf_init(n_a,2,N_lifetime)
double precision :: pmf(n_a,2,N_lifetime)
double precision :: L_supply
double precision :: K_supply
integer						:: 				pf_i_apr(n_a,2,N_lifetime)
!double precision 						:: pmf_init_flat(n_a,2,N_lifetime)
integer 										::			  i_stat   ! Used for writing data after program has run.
integer :: i_age 
integer :: i_exp
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

double precision :: mass(N_lifetime)
double precision :: t_mass = 0


! Discretizing the state space (capital)
do i_a = 1,n_a
	grid_a(i_a) = min_a + (dble(i_a) - 1d0)*step_a
end do

! Setting up Policy Function guesses
do i_a = 1,n_a
do i_age = 1,N_lifetime
	pf_c(i_a,1,i_age) 			    = 0d0
	pf_a(i_a,1,i_age) 		    	= 0d0
	pf_v(i_a,1,i_age)    	      = 0d0
	pf_l(i_a,1,i_age)    	      = 0d0
	pf_c(i_a,2,i_age) 			    = 0d0
	pf_a(i_a,2,i_age) 		    	= 0d0
	pf_v(i_a,2,i_age)    	      = 0d0
	pf_l(i_a,2,i_age)    	      = 0d0

	pmf_init(i_a,1,i_age) = 0d0
	pmf_init(i_a,2,i_age) = 0d0
	pf_i_apr(i_a,1,i_age) = 0
	pf_i_apr(i_a,2,i_age) = 0
end do
end do
! calc working mass
mass(1) = 1
t_mass = 1
do i_age = 2,N_lifetime
mass(i_age) = mass(i_age - 1)/(1+n_popg)
t_mass = t_mass + mass(i_age)
end do
!write(*,*) "mass: ", t_mass
working_mass = 0
do i_age = 1,J_retire-1
working_mass = working_mass+mass(i_age)/t_mass
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
double precision 						:: 				pf_c_up(n_a,n_z,N_lifetime)
double precision 						:: 				pf_a_up(n_a,n_z,N_lifetime)
double precision 						:: 				pf_v_up(n_a,n_z,N_lifetime)
double precision 						:: 				pf_l_up(n_a,n_z,N_lifetime)

double precision 						:: 				diff_c
double precision 						:: 				diff_a
double precision 						:: 				diff_v
double precision 						:: 				diff_l
double precision 						:: 				pmf_diff
double precision 						:: 				max_diff
double precision :: del
double precision :: l_today_temp
double precision :: l_today_max
double precision :: eff_w
double precision :: Ph
double precision :: Pl
double precision :: sum_pmf
double precision :: L_demand
double precision :: K_demand
double precision :: Z
integer :: apr_floor
integer :: next_e
integer :: next_u
integer :: decl
!integer :: i_exp
!integer :: n_samples = 10 ! number of samples for iterated grid search
integer :: i_sample
integer :: i_min(1)
integer :: i_mid 
integer :: i_n
integer :: n
integer :: n_mid = (n_samples+1)/2
integer :: age 
integer :: converged_outer = 0
integer :: i_apr_max
integer :: i_go_here
double precision :: apr_sum
double precision :: mu
double precision :: sig
double precision :: l_today
double precision :: mass_newborn
double precision :: popweights_true(N_lifetime)
double precision :: del_L
integer :: i_con
double precision :: kdiff
double precision :: pmf_prime(n_a,2)
double precision :: q_min = cBET!0.5!cBET!0.99!0.99!cBET - 0.1!0.4d0
double precision :: q_max = 1.0!2.0!0.5!1.0!cBET + 0.05!cBET
double precision :: q_grid(n_samples-2)
double precision :: apr_mid = 1d0
double precision :: age_prod( J_retire-1)
double precision :: L_init = 0.387!0.4
double precision :: K_init = 5.68!8 ! these are just guesses for now
double precision :: tune_KK 
double precision :: tune_LL

	  age_prod(1) = 0.59923239  ! I do not know how to read data into fortran yet lol
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


do i_exp = 1,2 ! First thing to do: calculate the steady states with the two different theta values
converged_outer = 0
if (i_exp<1.5) then
cTHETA = 0.11d0
cGAMMA = 0.42d0!1d0 ! labor is exogenous here
cZh = 3.0d0
L_init = 0.343207!0.7542561 
K_init =  3.35963!4.55470!3.64!7.2388
tune_K = 0.95 ! I know where these are so I am going to converge fast
tune_L = 0.95 ! what is going to happen is on the first run I will calculate L and then
			 ! it will stay there forever
else if (i_exp<2.5) then
cTHETA = 0.0d0
cGAMMA = 0.42d0!1d0 ! labor is exogenous here
cZh = 3.0d0
L_init = 0.36532!0.7542561 
K_init =  4.60341!9.89349
tune_K = 0.95
tune_L = 0.95 
end if

tune_KK = 0.5d0
tune_LL = 0.5d0
! write  loop out here and try to compute
L_demand = L_init
K_demand = K_init

do while (converged_outer==0)
rental = cALPHA* ((K_demand)**(cALPHA - 1))*((L_demand)**(1 - cALPHA)) - cDELTA !note: delta here!
wage = (1-cALPHA)* ((K_demand)**(cALPHA))*((L_demand)**(-1*cALPHA))
benefits = cTHETA*wage*L_demand/(1-working_mass) ! none of this changes

converged = 0
it = 1
! We use backwards induction here.
do i_age = 1,N_lifetime 
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
		if (age == N_lifetime) then ! we eat our budget
			c_today = y_today
			if (c_today>0d0) then
			v_today = (c_today ** ((1-cSIGMA)*cGAMMA))/(1-cSIGMA)
			else
			v_today = -1d10
			end if
			a_tomorrow = 0
			l_today_max = 0d0
		else if (age>=J_retire) then
			i_apr = apr_floor
			decl = 0
			do while (decl == 0 .and. i_apr<=n_a)
				a_tomorrow = grid_a(i_apr)
				c_today_temp = y_today - a_tomorrow
				if (c_today_temp>0d0) then 
					if (i_z < 1.5) then
					v_tomorrow = cPhh*pf_v(i_apr,1,age+1) + cPhl*pf_v(i_apr,2,age+1) 
				else
					v_tomorrow = cPlh*pf_v(i_apr,1,age+1) + cPll*pf_v(i_apr,2,age+1) 
				end if
					v_today_temp = (c_today_temp**((1-cSIGMA)*cGAMMA))/(1-cSIGMA) + cBET*v_tomorrow
					if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
						v_today = v_today_temp
						a_tomorrow_max = a_tomorrow
						apr_floor = i_apr
						pf_i_apr(i_a,i_z,age) = i_apr
					else
						decl = 1
					end if
				end if
				i_apr = i_apr+1
			end do
			a_tomorrow = a_tomorrow_max
			c_today = y_today - a_tomorrow
			l_today_max = 0d0
		else	!working years
		i_apr = apr_floor
		decl = 0 	
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
			c_today_temp = eff_w*l_today_temp+ (1+rental)*a_today - a_tomorrow
			if (c_today_temp>0d0) then 	! if this thing is negative then don't check it
			if (i_z <1.5) then
				v_tomorrow = cPhh*pf_v(i_apr,1,age+1) + cPhl*pf_v(i_apr,2,age+1)  							
			else
				v_tomorrow = cPlh*pf_v(i_apr,1,age+1) + cPll*pf_v(i_apr,2,age+1)  
			end if
			v_today_temp = ((c_today_temp)**(cGAMMA)*(1-l_today_temp)**(1-cGAMMA) )**(1-cSIGMA)/(1-cSIGMA) + cBET*v_tomorrow

			if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
				v_today = v_today_temp
				a_tomorrow_max = a_tomorrow
				l_today_max = l_today_temp
				apr_floor = i_apr
				pf_i_apr(i_a,i_z,age) = i_apr
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
	  	pf_c_up(i_a,i_z,age) = c_today 												
	  	pf_a_up(i_a,i_z,age) = a_tomorrow												
	  	pf_v_up(i_a,i_z,age) = v_today	
	  	pf_l_up(i_a,i_z,age) = l_today_max									
	end do
	end do
	pf_c 		= pf_c_up
	pf_a 		= pf_a_up
	pf_v			= pf_v_up
	pf_l			= pf_l_up
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

sum_pmf = 0d0 ! need to make pmf sum to 1
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
mass_newborn = pmf(1,1,1)+pmf(1,2,1)

! need to calculate labor supply, capital supply
L_supply = 0
do i_z = 1,2
if (i_z<1.5) then
Z = cZh
else
Z = cZl
end if
do i_age = 1,(J_retire - 1)
	do i_a = 1,n_a
		if (pmf(i_a,i_z,i_age)>0d0) then
			L_supply = L_supply + pmf(i_a,i_z,i_age)*Z*age_prod(i_age)*pf_l(i_a,i_z,i_age)
		end if
	end do
end do
end do
K_supply = 0
do i_z = 1,2
do i_age = 1,N_lifetime
	do i_a = 1,n_a
	if (pmf(i_a,i_z,i_age)>0d0) then
		K_supply = K_supply + pmf(i_a,i_z,i_age)*grid_a(i_a)!pmf(i_a,i_z,i_age)*pf_a(i_a,i_z,i_age)
	end if
	end do
end do
end do

diff_a = abs(K_supply - K_demand)
diff_l = abs(L_supply - L_demand)
if ((diff_a<a_tol).and.(diff_l<l_tol)) then
converged_outer = 1
else
K_demand = tune_K *K_demand + (1-tune_K) * K_supply
L_demand = tune_L *L_demand + (1-tune_L) * L_supply
!write(*,*) "K supply: ", K_supply
!write(*,*) "K demand: ", K_demand
!write(*,*) "L supply: ", L_supply
!write(*,*) "L demand: ", L_demand
end if
end do


if (i_exp<1.5) then
!if (i_exp>1.5) then
K_0ss = K_supply
L_0ss = L_supply
do i_a = 1,n_a
do i_z = 1,n_z
do i_age = 1,N_lifetime
tran_pmf(i_a,i_z,i_age,1) = pmf(i_a,i_z,i_age)
end do
end do
end do
do n = 1,N_transition
do i_z = 1,n_z
tran_pmf(1,i_z,1,n) = pmf(1,i_z,1) !newborn mass in all periods 
end do
!write(*,*) "Newborn mass: ", (tran_pmf(1,1,1,n)+tran_pmf(1,2,1,n)) ! all correct here
end do

do i_age = 1,N_lifetime
sum_pmf = 0d0
do i_a = 1,n_a
do i_z = 1,n_z
sum_pmf = sum_pmf + pmf(i_a,i_z,i_age)
end do
end do
popweights_true(i_age) = sum_pmf
!if (i_age>1.5) then
!write(*,*) "implied growth rate: ", popweights_true(i_age)/popweights_true(i_age-1)
!end if
end do

else
K_Tss = K_supply
L_TSS = L_supply
end if
end do
write(*,*) "K_0ss: ", K_0ss
write(*,*) "L_0ss: ", L_0ss
write(*,*) "K_Tss: ", K_Tss
write(*,*) "L_Tss: ", L_Tss

! debugging 
!K_0ss = K_Tss
!L_0ss = L_Tss

do i_a = 1,n_a
do i_z = 1,n_z
do i_age = 1,N_lifetime
tran_pf_v(i_a,i_z,i_age,N_transition+1) = pf_v(i_a,i_z,i_age) ! final tran_pf_v is new ss
end do
end do
end do

do i_exp = 1,3!,3
! initialize guess of transition to capital
del = (K_Tss - K_0ss)/N_transition
del_L = (L_Tss - L_0ss)/N_transition
do i_a = 1,N_transition
	K_0transition(i_a) = K_0ss + (i_a-1) * del ! first k_0transition is at K_0ss because at the moment of MIT shock the households are in old SS asset positions
	L_0transition(i_a) = L_0ss + (i_a-1) * del_L
end do

! ----------------------------------------------------------
! HERE I DO THE TRANSITION PATH STUFF
! ----------------------------------------------------------
!! NEW LOOP - BACKWARDS INDUCT AT EACH STAGE TO GET EQM


converged_outer = 0
do while (converged_outer == 0)
n = N_transition + 1 - i_n
do i_n = 1,N_transition
if (i_exp<2.5) then
cTHETA = 0d0
else
if (n<22) then !dean's indexing starts at 0, mine starts at 1
cTHETA = 0.11d0
else
cTHETA = 0d0
end if
end if
K_demand = K_0transition(n)
L_demand = L_0transition(n)
! calc prices
rental = cALPHA* ((K_demand)**(cALPHA - 1))*((L_demand)**(1 - cALPHA)) - cDELTA !note: delta here!
wage = (1-cALPHA)* ((K_demand)**(cALPHA))*((L_demand)**(-1*cALPHA))
benefits = cTHETA*wage*L_demand/(1-working_mass) ! none of this changes
! iterate backwards to get transition rule
do i_age = 1,N_lifetime
age = N_lifetime + 1 - i_age
do i_z = 1,2
v_today = -10d10
if (i_z<1.5) then
Z = cZh
else
Z = cZl
end if
apr_floor = 1	
do i_a = 1,n_a
a_today = grid_a(i_a)
if (age == N_lifetime) then ! old, about to die
c_today = (1+rental)* a_today + benefits
if (c_today>1d-10) then
v_today = (c_today**((1-cSIGMA)*cGAMMA))/(1- cSIGMA)
else
v_today = -10d10
end if
l_today = 0
a_tomorrow_max = 0
i_apr_max = 0
else if (age >= J_retire) then ! old, will live another period
			i_apr = apr_floor
			decl = 0
			do while (decl == 0 .and. i_apr<=n_a)
				a_tomorrow = grid_a(i_apr)
				c_today_temp = (1+rental)*a_today + benefits - a_tomorrow
				if (c_today_temp>0d0) then 
					if (i_z < 1.5) then
					v_tomorrow = cPhh*tran_pf_v(i_apr,1,age+1,n+1) + cPhl*tran_pf_v(i_apr,2,age+1,n+1) ! <-- this needs to change tran_pf_v(i_apr,1,age+1,n+1)
					else
					v_tomorrow = cPlh*tran_pf_v(i_apr,1,age+1,n+1) + cPll*tran_pf_v(i_apr,2,age+1,n+1) ! <-- this needs to change
					end if
					v_today_temp = (c_today_temp**((1-cSIGMA)*cGAMMA))/(1-cSIGMA) + cBET*v_tomorrow
					if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
						v_today = v_today_temp
						a_tomorrow_max = a_tomorrow
						apr_floor = i_apr
						!pf_i_apr(i_a,i_z,age) = i_apr
						i_apr_max = i_apr
					else
						decl = 1
					end if
				end if
				i_apr = i_apr+1
			end do
			a_tomorrow = a_tomorrow_max
			c_today = (1+rental)*a_today + benefits - a_tomorrow
			l_today_max = 0d0

else !! working age 
! loop through a' possibilities and pick

		i_apr = apr_floor
		decl = 0 	
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
				v_tomorrow = cPhh*tran_pf_v(i_apr,1,age+1,n+1) + cPhl*tran_pf_v(i_apr,2,age+1,n+1)					
			else
				v_tomorrow = cPlh*tran_pf_v(i_apr,1,age+1,n+1) + cPll*tran_pf_v(i_apr,2,age+1,n+1)  
			end if
			v_today_temp = ((c_today_temp)**(cGAMMA)*(1-l_today_temp)**(1-cGAMMA) )**(1-cSIGMA)/(1-cSIGMA) + cBET*v_tomorrow

			if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
				v_today = v_today_temp
				a_tomorrow_max = a_tomorrow
				l_today_max = l_today_temp
				apr_floor = i_apr
				!pf_i_apr(i_a,i_z,age) = i_apr
				i_apr_max = i_apr
			else
				decl = 1
			end if
			end if
			i_apr = i_apr+1
		end do
		a_tomorrow = a_tomorrow_max
		c_today = eff_w*l_today_max + (1+rental)*a_today - a_tomorrow

end if
! okay now I have my choice variables. 

tran_pf_v(i_a,i_z,age,n) = v_today
tran_pf_a(i_a,i_z,age,n) = a_tomorrow_max
tran_pf_c(i_a,i_z,age,n) = c_today
tran_pf_i_apr(i_a,i_z,age,n) = i_apr_max
tran_pf_l(i_a,i_z,age,n) = l_today_max
end do
end do
end do
end do
! okay now I have transition functions, now I need to figure out where ppl are

! reset relevant parts of tran pmf to zero before we iterate forwards
do i_n = 2,N_transition ! don't change first transition pmf because people have steady state assets at that time
do i_a = 1,n_a
do i_z = 1,n_z
do i_age = 2,N_lifetime
tran_pmf(i_a,i_z,i_age,i_n) = 0d0
end do
end do
end do
end do

do i_n = 1,N_transition !newborns never have mass outside of zero asset position
do i_a = 2,n_a
do i_z = 1,n_z
tran_pmf(i_a,i_z,1,i_n) = 0d0
end do
end do
end do

!do i_n = 1,N_transition 
!sum_pmf = 0d0
!do i_a = 1,n_a ! for each current asset holdings ...
!do i_z = 1,n_z ! for each productivity state ...
!do i_age = 1,1!1,(N_lifetime) ! for each current age that will be alive next period ...
!sum_pmf = sum_pmf + tran_pmf(i_a,i_z,i_age,i_n)
!if (i_a>1 .and. tran_pmf(i_a,i_z,i_age,i_n)>0d0) then
!write(*,*) "UNEXPECTED MASS AT GRID POINT i_a: ", i_a
!write(*,*) "UNEXPECTED MASS AT GRID POINT i_z: ", i_z
!end if 
!end do
!end do
!end do
!write(*,*) "mass sums to ", sum_pmf
!if (abs(sum_pmf - 1)>1d-9) then
!write(*,*) "--------------------------"
!write(*,*) "time ", i_n
!write(*,*) "young mass sums to ", sum_pmf
!write(*,*) "--------------------------"
!end if
!end do


! iterate forwards to get capital levels
do i_n = 1,(N_transition-1) ! for each point in the transition path ...
do i_a = 1,n_a ! for each current asset holdings ...
do i_z = 1,n_z ! for each productivity state ...
do i_age = 1,(N_lifetime-1) ! for each current age that will be alive next period ... (note: N_lifetime will be dead next period so their mass drops out)
! figure out where they will be next period.
if (tran_pmf(i_a,i_z,i_age,i_n)>0d0) then ! do we have any mass to move forwards? If not, nobody exists here so position is irrelevant
i_go_here = tran_pf_i_apr(i_a,i_z,i_age,i_n)
if (i_z<1.5) then
tran_pmf(i_go_here,1,i_age+1,i_n+1) = tran_pmf(i_go_here,1,i_age+1,i_n+1) &
															+ tran_pmf(i_a,i_z,i_age,i_n) * cPhh * (1/(1+n_popg)) 		! mass I have * decay * prob good state|s
tran_pmf(i_go_here,2,i_age+1,i_n+1) = tran_pmf(i_go_here,2,i_age+1,i_n+1)  &
															+ tran_pmf(i_a,i_z,i_age,i_n) * cPhl * (1/(1+n_popg))
else
tran_pmf(i_go_here,1,i_age+1,i_n+1) = tran_pmf(i_go_here,1,i_age+1,i_n+1) &
															+ tran_pmf(i_a,i_z,i_age,i_n) * cPlh * (1/(1+n_popg))
tran_pmf(i_go_here,2,i_age+1,i_n+1) = tran_pmf(i_go_here,2,i_age+1,i_n+1) &
															+ tran_pmf(i_a,i_z,i_age,i_n) * cPll * (1/(1+n_popg))
if ((tran_pf_i_apr(i_a,i_z,i_age,i_n)<1).and.i_age<65) then !never tripped
write(*,*) "warning: policy function incorrect at i_n ",i_n
write(*,*) "warning: policy function incorrect at i_a ",i_a
write(*,*) "warning: policy function incorrect at i_z ",i_z
write(*,*) "warning: policy function incorrect at i_age ",i_age
end if
end if
end if
end do
end do
end do
! Renormalize (machine precision issues can cause pmf to not sum to one)
do i_age = 1,(N_lifetime) ! for each current age that will be alive next period ... 
sum_pmf = 0d0
do i_a = 1,n_a ! for each current asset holdings ...
do i_z = 1,n_z ! for each productivity state ...
sum_pmf = sum_pmf + tran_pmf(i_a,i_z,i_age,i_n+1)
end do
end do
do i_a = 1,n_a ! for each current asset holdings ...
do i_z = 1,n_z ! for each productivity state ...
 tran_pmf(i_a,i_z,i_age,i_n+1) = popweights_true(i_age)*tran_pmf(i_a,i_z,i_age,i_n+1)/sum_pmf
end do
end do
end do
!do i_a = 1,n_a ! for each current asset holdings ...
!do i_z = 1,n_z ! for each productivity state ...
!do i_age = 1,(N_lifetime) ! for each current age that will be alive next period ... (note: N_lifetime will be dead next period so their mass drops out)
!tran_pmf(i_a,i_z,i_age,i_n+1) = tran_pmf(i_a,i_z,i_age,i_n+1)/sum_pmf
!end do
!end do
!end do
end do


do i_n = 1,N_transition ! for each point on the transition path ...
K_1transition(i_n) = 0d0
L_1transition(i_n) = 0d0
! compute capital supply
do i_a = 1,n_a
do i_z = 1,n_z
if (i_z < 1.5) then
Z = cZh
else
Z = cZl
end if
do i_age = 1,N_lifetime
K_1transition(i_n) = K_1transition(i_n) + tran_pmf(i_a,i_z,i_age,i_n)*grid_a(i_a)
L_1transition(i_n) = L_1transition(i_n) + tran_pmf(i_a,i_z,i_age,i_n)*age_prod(i_age)*Z*tran_pf_l(i_a,i_z,i_age,i_n) ! pmf(i_a,i_z,i_age)*Z*age_prod(i_age)*pf_l(i_a,i_z,i_age)
end do
end do
end do
end do

! check convergence of capital supply with capital demand
diff_a = maxval(abs(K_1transition - K_0transition))
diff_l = maxval(abs(L_1transition - L_0transition))
if ((diff_a<tol).and.(diff_l<tol)) then
converged_outer = 1
else
write(*,*) "Max diff: ", diff_a
do i_n = 1,N_transition
K_0transition(i_n) = tune_KK*K_0transition(i_n) + (1-tune_KK)*K_1transition(i_n)
L_0transition(i_n) = tune_LL*L_0transition(i_n) + (1-tune_LL)*L_1transition(i_n)
!converged_outer = 1
end do
end if
end do
! display how close we got to final tolerance level

write(*,*) "Final K diff: ", K_1transition(N_transition)-K_Tss
write(*,*) "Final L diff: ", L_1transition(N_transition)-L_Tss
if (i_exp<1.5) then
! figure out "new steady state"
kdiff = 10d0
i_con = 0
do i_n = 1,(N_transition-1)
if (abs(K_1transition(i_n+1)-K_1transition(i_n))<kdiff) then
kdiff = abs(K_1transition(i_n+1)-K_1transition(i_n))
i_con = i_n
end if
end do


do i_a = 1,n_a
do i_z = 1,n_z
do i_age = 1,N_lifetime
tran_pf_v(i_a,i_z,i_age,N_transition+1) = tran_pf_v(i_a,i_z,i_age,i_con)
end do
end do
end do

K_Tss = K_1transition(i_con)
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
open(unit = 2, file = 'pfs_K.dat', status = 'replace', action = 'write', iostat = i_stat)
200 format(f25.15,2x,f25.15,2x,f25.15,2x) !,f25.15,2x,f25.15,2x

do i_exp = 1,N_transition
     write(2,200) K_0transition(i_exp), K_0ss, K_Tss !res_K(i_exp), res_L(i_exp), res_w(i_exp), res_r(i_exp), res_b(i_exp), res_welf(i_exp), res_cv(i_exp) !pf_v(i_a,50), pf_a(i_a,20), pf_a_b(i_a,20), pf_a(i_a,1), pf_a_b(i_a,2), &
      !pmf(i_a,1,N_lifetime), pmf(i_a,2,N_lifetime)
     !, pf_i_apr(i_a,1), pf_i_apr(i_a,2)!grid_a(i_a), pf_c(i_a), pf_a(i_a), pf_v(i_a), pmf(i_a,1), pmf(i_a,2), pmf_init_flat(i_a,1), pmf_init_flat(i_a,2)!, pf_i_apr(i_a,1), pf_i_apr(i_a,2)
end do

return

end subroutine coda
