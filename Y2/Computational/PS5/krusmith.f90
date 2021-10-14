! ************************************************************************
! 
! Author : Michael Nattinger
!
! Implements Krusell and Smith (1998)
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
integer, parameter ::T_timeseries = 11000
integer, parameter :: T_burn = 1000 ! burn-in
integer, parameter :: T_net = T_timeseries - T_burn
integer, parameter :: N_shocks = 50 ! retirement age
double precision :: shock_grid(T_timeseries,N_shocks) 
double precision :: a0 = 0.095d0
double precision :: a1 = 0.9d0!0.999d0 !initial guesses
double precision :: b0 = 0.085d0
double precision :: b1 = 0.9d0!0.999d0
double precision, parameter :: k_lower = 0d0
double precision, parameter :: KK_lower = 11d0
double precision, parameter :: k_upper = 15d0
double precision, parameter :: KK_upper = 15d0!15d0
integer, parameter :: n_k = 100
integer, parameter :: n_KK = 100 ! 100 each
double precision :: grid_k(n_k)
double precision :: grid_KK(n_KK)
double precision :: step_k =(k_upper - k_lower)/(n_k - 1)
double precision :: step_KK =(KK_upper - KK_lower)/(n_KK - 1)
integer, parameter :: n_z = 2
integer, parameter :: n_e = 2 ! number of employment statuses
double precision :: KK_timeseries(T_timeseries)
double precision :: k_timeseries(T_timeseries,N_shocks)
integer :: z_timeseries(T_timeseries)
integer :: e_timeseries(T_timeseries,N_shocks)

double precision :: R2
double precision :: OLS(4,1) 
double precision, parameter :: res_tol = 2d-2 ! tolerance for OLS coefficients
double precision :: OLS_lag(4,1) ! tolerance for OLS coefficients
double precision, parameter :: K_SS = 5.7163
double precision, parameter :: cBET 				= 0.99d0 	! Discount Factor
double precision, parameter :: cALPHA 			= 0.36d0	! Capital Share
double precision, parameter :: cDELTA = 0.025
!double precision, parameter :: cERGO_g = 0.2037
!double precision, parameter :: cERGO_b = 1-cERGO_g ! ergodic distribution
double precision, parameter :: cZh 				= 1.01d0	! High productivity
double precision, parameter :: cZl  				= 0.99d0	! Low productivity
double precision, parameter :: cEPSh 				= 0.3271d0	! High empl
double precision, parameter :: cEPSl  				= 0d0	! Low empl
!double precision, parameter :: cPhh 				= 0.9261d0	! Probability of good -> good transition
!double precision, parameter :: cPhl 				= 1-cPhh	! Probability of good -> bad transition
!double precision, parameter :: cPll 				= 0.9811d0	! Probability of bad -> good transition
!double precision, parameter :: cPlh 				= 1-cPll	! Probability of bad -> bad transition
! still define these here but not as a parameter - these will be edited in the loop as appropriate
double precision :: wage 
double precision :: rental 
double precision :: benefits
double precision :: Pr(4,4)
double precision :: uh,ul

double precision :: pgg
double precision :: pbb
double precision :: pbg
double precision :: pgb
double precision :: Prgg
double precision :: Prbb
double precision :: Prbg
double precision :: Prgb
double precision :: diff


! Tolerance level for convergence and max itations
integer, parameter 	:: 	max_it 		= 10000 ! Maximum Number of Iterations
integer :: it = 1 		! Itation Counter
integer :: converged		= 0			! Dummy for VFI Convergence

! -----------------------------------------------------------------------
! ****************************GRID SET**********************************
! -----------------------------------------------------------------------
! Set up for discritizing the state space 
integer				  :: 				i_a, i_apr		
integer :: i_k, i_KK							! Ieration Counters for k_today and k_tomorrow grid
!integer, parameter 				  :: 				n_a 				= 1000!500!2000																! Size of k grid
!integer, parameter :: n_z = 2
!double precision 						:: 				grid_a(n_a)																				! Allocate Space for k grid
!double precision, parameter :: 				min_a 			= 0d0															! Minimum of k grid
!double precision, parameter :: 				max_a 			= 40d0																! Maximum of k grid
!double precision, parameter :: 				step_a 			= (max_a - min_a)/(dble(n_a) - 1d0) 	! Step of k grid
double precision					  :: 				a_today
double precision					  :: 				a_tomorrow
integer :: i_z ! 1 good state, 2 bad state
integer :: n_good

! Global variables for Dynamic Progamming
double precision 						:: 				c_today
double precision 						:: 				c_today_temp
double precision 						:: 				y_today
double precision 						:: 				a_tomorrow_max
double precision 						:: 				v_today
double precision 						:: 				v_today_temp
double precision 						:: 				v_tomorrow


! Allocating space for Policy Functions
double precision 						:: 				pf_c(n_k,n_z,n_KK,n_e) ! Allocate Space for Conumption Policy Function ! no _b means it is implicitly _g
double precision 						:: 				pf_a(n_k,n_z,n_KK,n_e) ! Allocate Space for Capital Policy Function
double precision 						:: 				pf_v(n_k,n_z,n_KK,n_e) ! Allocate Space for Value Function
double precision 						:: 				pf_l(n_k,n_z,n_KK,n_e) ! Allocate Space for Value Function

double precision :: L_supply
double precision :: K_supply
integer						:: 				pf_i_apr(n_k,n_z,n_KK,n_e)
!double precision 						:: pmf_init_flat(n_a,2,N_lifetime)
integer 										::			  i_stat   ! Used for writing data after program has run.
integer :: i_e 
end module params_grid


! ************************************************************************
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------
! program : KruSmith
!
! Description : Implements Krusell and Smith (1998) - Michael Nattinger
! ------------------------------------------------------------------------
! ------------------------------------------------------------------------

program KruSmith

use params_grid


implicit none

! Begin Computational Timer
integer 										::				beginning, end, rate

call system_clock(beginning, rate)

! Initialize Grids and Policy Functions
call housekeeping()

! Do Value Function Iteration & Search for Fixed Point
!call bellman(n_good)
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

end program KruSmith

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

double precision :: durug ! from 
double precision :: unempg
double precision :: durgd
double precision :: unempb
double precision :: durbd
double precision :: durub
double precision :: pgg00
double precision :: pbb00
double precision :: pbg00
double precision :: pgb00
double precision :: pgg01
double precision :: pbb01
double precision :: pbg01
double precision :: pgb01
double precision :: pgg10
double precision :: pbb10
double precision :: pbg10
double precision :: pgb10
double precision :: pgg11
double precision :: pbb11
double precision :: pbg11
double precision :: pgb11
double precision :: rn
!double precision :: Prgg
!double precision :: Prgb
!double precision :: Prbg
!double precision :: Prbb
integer :: col
integer :: i_hh
integer :: i_sim
integer :: row_1
integer :: row_2

! Discretizing the state space (capital)
do i_k = 1,n_k
	grid_k(i_k) = k_lower + (dble(i_k) - 1d0)*step_k
end do

do i_KK = 1,n_KK
	grid_KK(i_KK) = KK_lower + (dble(i_KK) - 1d0)*step_KK
end do


! Setting up Policy Function guesses
do i_k = 1,n_k
do i_KK = 1,n_KK
do i_z = 1,2
do i_e = 1,2
	pf_c(i_k,i_z,i_KK,i_e) 			    = 0d0
	pf_a(i_k,i_z,i_KK,i_e)  		    	= 0d0
	pf_v(i_k,i_z,i_KK,i_e)    	      = 0d0
end do
end do
end do
end do


!%parameters of transition matrix:
durug=1.5
unempg=0.04
durgd=8.0
unempb=0.1
durbd=8.0
durub=2.5
!%transition probabilities
pgg00 = (durug-1)/durug
pbb00 = (durub-1)/durub
pbg00 = 1.25*pbb00
pgb00 = 0.75*pgg00
pgg01 = (unempg - unempg*pgg00)/(1-unempg)
pbb01 = (unempb - unempb*pbb00)/(1-unempb)
pbg01 = (unempb - unempg*pbg00)/(1-unempg)
pgb01 = (unempg - unempb*pgb00)/(1-unempb)
pgg = (durgd-1)/durgd
pgb = 1 - (durbd-1)/durbd
pgg10 = 1 - (durug-1)/durug
pbb10 = 1 - (durub-1)/durub
pbg10 = 1 - 1.25*pbb00
pgb10 = 1 - 0.75*pgg00
pgg11 = 1 - (unempg - unempg*pgg00)/(1-unempg)
pbb11 = 1 - (unempb - unempb*pbb00)/(1-unempb)
pbg11 = 1 - (unempb - unempg*pbg00)/(1-unempg)
pgb11 = 1 - (unempg - unempb*pgb00)/(1-unempb)
pbg = 1 - (durgd-1)/durgd
pbb = (durbd-1)/durbd
!%matrix
pr(1,1) = pgg*pgg11
pr(2,1) = pbg*pbg11
pr(3,1) = pgg*pgg01
pr(4,1) = pbg*pbg01
pr(1,2) = pgb*pgb11
pr(2,2) = pbb*pbb11
pr(3,2) = pgb*pgb01
pr(4,2) = pbb*pbb01
pr(1,3) = pgg*pgg10
pr(2,3) = pbg*pbg10
pr(3,3) = pgg*pgg00
pr(4,3) = pbg*pbg00
pr(1,4) = pgb*pgb10
pr(2,4) = pbb*pbb10
pr(3,4) = pgb*pgb00
pr(4,4) = pbb*pbb00
Prgg = pr(1,1) + pr(1,3)
Prgb = 1 - Prgg
Prbb = pr(4,4) + pr(4,2)
Prbg = 1-Prbb

! when I initialize shocks don't forget to count n_g

! generate productivity time series
z_timeseries(1) = 1
n_good = 0
do i_sim = 2,T_timeseries
call random_number(rn)
! I need to determine markov transition probabilities 
if (z_timeseries(i_sim-1) < 1.5) then ! good state
if (rn>Prgg) then ! transition to bad
z_timeseries(i_sim) = 2
else
z_timeseries(i_sim) = 1
if (i_sim>T_burn) then
n_good = n_good + 1
end if
end if
else
if (rn>Prbg) then ! transition to bad
z_timeseries(i_sim) = 2
else
z_timeseries(i_sim) = 1
if (i_sim>T_burn) then
n_good = n_good + 1
end if
end if
end if
end do

do i_hh = 1,N_shocks ! generate employment status time series
! initial position
call random_number(rn)
if (rn> uh) then !employed
e_timeseries(1,i_hh) = 1
else ! unemployed
e_timeseries(1,i_hh) = 2
end if
! finish
do i_sim = 2,T_timeseries
call random_number(rn) ! get my uniform draw
if (z_timeseries(i_sim - 1)<1.5) then ! I was in good state
if (e_timeseries(i_sim - 1, i_hh)<1.5) then ! I was employed
col = 1
else ! I was unemployed
col = 3
end if
else ! I was in bad state
if (e_timeseries(i_sim - 1, i_hh)<1.5) then ! I was employed
col = 2
else ! I was unemployed
col = 4
end if
end if
if (z_timeseries(i_sim )<1.5) then
row_1 = 1
row_2 = 3
else
row_1 = 2
row_2 = 4
end if
if (rn<(Pr(row_1,col)/(Pr(row_1,col)+Pr(row_2,col)))) then ! I am employed P(e'=1|e,z,z')
e_timeseries(i_sim,i_hh) = 1
else
e_timeseries(i_sim,i_hh) = 2
end if
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

!subroutine bellman(n_g)
subroutine bellman()
use params_grid


implicit none
!integer, intent(in) :: n_g
!integer, parameter :: n_g = n_good
!integer, parameter :: n_b = T_timeseries-T_burn - n_g
!double precision :: good_resid(n_g,1)
!double precision :: good_Y(n_g,1)
!double precision :: good_X(n_g,2)
!double precision :: good_OLS(2,1)
!double precision :: bad_Y(n_b,1)
!double precision :: bad_resid(n_b,1)
!double precision :: bad_X(n_b,2)
!double precision :: bad_OLS(2,1)
double precision :: YY(T_net,1)
double precision :: XX(T_net,4)
double precision :: resid(T_net,1)
integer :: decl, i_g, i_hh, i_b, i_kkpr, i_sim
double precision :: k_today, KK_tomorrow

! allocating space for policy function updates
double precision 						:: 				pf_c_up(n_k,n_z,n_KK,n_e)
double precision 						:: 				pf_a_up(n_k,n_z,n_KK,n_e)
double precision 						:: 				pf_v_up(n_k,n_z,n_KK,n_e)
double precision 						:: 				pf_l_up(n_k,n_z,n_KK,n_e)
double precision :: iXX(4,4)
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
double precision :: K_demand, meanlogkk
double precision :: V_kk(n_KK)
double precision :: Z, MSE, Rsq, SSE, KK_expected
integer :: apr_floor
integer :: next_e
integer :: next_u
!integer :: decl
!integer :: i_exp
!integer :: n_samples = 10 ! number of samples for iterated grid search
integer :: i_sample
integer :: i_min(1)
integer :: i_mid 
integer :: i_n,n_b, e
integer :: n
integer :: age 
integer :: converged_outer = 0
integer :: i_apr_max
integer :: colnum

converged_outer = 0 ! step 1: value function iteration | ols parameters
OLS_lag(2,1) = a1
OLS_lag(1,1) = a0
OLS_lag(4,1) = b1
OLS_lag(3,1) = b0

uh = 0.04d0
ul = 0.1d0

do while (converged_outer==0)
converged = 0
it = 1
do while(converged==0)
	!write(*,*) "beginning iteration on VFI: ",it
	do i_KK = 1,n_KK
	K_supply = grid_KK(i_KK)
	do i_z = 1,2 ! Aggregate productivity
	if (i_z<1.5) then
		L_supply = cEPSh*(1-uh)
		Z = cZh
	else
		L_supply = cEPSh*(1-ul)
		Z = cZl
	end if
	wage = (1-cALPHA)*Z*(K_supply/L_supply)**(cALPHA)
	rental = (cALPHA)*Z*(K_supply/L_supply)**(cALPHA-1)
	!write(*,*) "rental", rental
	if (i_z<1.5) then
		KK_tomorrow = exp(a0 + a1*log(K_supply))
	else
		KK_tomorrow = exp(b0 + b1*log(K_supply))
	end if
	do i_e = 1,2 ! employment level
	if (i_e<1.5) then
	if (i_z<1.5) then
	colnum = 1
	else
	colnum = 2
	end if
	else
	if (i_z<1.5) then
	colnum = 3
	else
	colnum = 4
	end if
	end if
	e = 2-i_e 
	apr_floor = 1						
	do i_k = 1,n_k ! Loop over all possible capital states	
		! ***********************************
		! Define the today's state variables
		! ***********************************
		k_today = grid_k(i_k)

		! ******************************************************
		! Solve for the optimal consumption / capital investment
		! ******************************************************
		v_today = -1d10 ! Set a very low bound for value
		i_apr = apr_floor
		decl = 0 	
		a_tomorrow_max =0d0
		y_today = (1 - cDELTA + rental)*k_today + wage*e*cEPSh
		do while (decl == 0 .and. i_apr<=n_k)				
			a_tomorrow = grid_k(i_apr)
			c_today_temp = y_today - a_tomorrow
			if (c_today_temp>0d0) then 	! if this thing is negative then don't check it
			do i_KKpr = 1,n_KK
				V_kk(i_KKpr) = Pr(1,colnum) *pf_v(i_apr,1,i_KKpr,1) & !<-- fill in first index from phil's response
							 + Pr(3,colnum) *pf_v(i_apr,1,i_KKpr,2) &
							 + Pr(2,colnum) *pf_v(i_apr,2,i_KKpr,1) &
							 + Pr(4,colnum) *pf_v(i_apr,2,i_KKpr,2)
							 if (abs(Pr(1,colnum)+Pr(2,colnum)+Pr(3,colnum)+Pr(4,colnum)-1)>1d-6) then
							 write(*,*) "sum: ", Pr(1,colnum)+Pr(2,colnum)+Pr(3,colnum)+Pr(4,colnum)
							 end if
							 !write(*,*) "sum: ", Pr(1,colnum)+Pr(2,colnum)+Pr(3,colnum)+Pr(4,colnum)
			end do
			call lin_interp_1d(grid_KK,KK_tomorrow,n_KK,V_KK,v_tomorrow) !calculate v_tomorrow
			!if (i_apr==85) then
			!write(*,*) "v_tomorrow",v_tomorrow
			!end if
			v_today_temp = log(c_today_temp) + cBET*v_tomorrow
			if (v_today_temp >= v_today) then ! if "temp" value is best so far, record value and capital choice
				v_today = v_today_temp
				a_tomorrow_max = a_tomorrow
				apr_floor = i_apr
				pf_i_apr(i_a,i_z,i_KK,i_e) = i_apr
			else
				decl = 1 ! we have begun declining and don't need to check any more a' values
			end if
			else
			decl = 1
			end if
			i_apr = i_apr+1
		end do
		a_tomorrow = a_tomorrow_max
		c_today = y_today - a_tomorrow
   		! *******************************
	  	! ****Update Policy Functions****
	  	! *******************************
	  	pf_c_up(i_k,i_z,i_KK,i_e) = c_today 												
	  	pf_a_up(i_k,i_z,i_KK,i_e) = a_tomorrow												
	  	pf_v_up(i_k,i_z,i_KK,i_e) = v_today	
	  	!if (a_tomorrow == k_upper) then
	  	!write(*,*) "a_tomorrow is at max:", k_upper
	  	!end if							
	end do
	end do
	end do
	end do
	! check for convergence
	diff = maxval(abs(pf_v - pf_v_up))
	pf_c = pf_c_up
	pf_a = pf_a_up
	pf_v = pf_v_up
	if (abs(diff)>res_tol) then
		!converged = 1
		!write(*,*) "not converged, abs(diff)", abs(diff)
	else
	converged = 1 ! just for writing code
	end if
	it = it+1
	!write(*,*) "abs(diff)", abs(diff)
end do
write(*,*) "abs(diff)", abs(diff)

! step 2: Simulate
write(*,*) "step"
do i_hh = 1,N_shocks ! initialize simulation capital levels
k_timeseries(1,i_hh) = K_SS
end do
do i_sim = 2,T_timeseries! for each simulation period...
if (z_timeseries(i_sim)<1.5) then
KK_expected = exp(a0+a1*log(KK_timeseries(i_sim-1)))
!write(*,*) "KK_expected", KK_expected
else
KK_expected = exp(b0+b1*log(KK_timeseries(i_sim-1)))
write(*,*) "KK_expected", KK_expected
end if
do i_hh = 1,N_shocks! for each household...
! use policy function to determine the policy of the household
call lin_interp_2d(grid_k,grid_KK,k_timeseries(1,i_hh),KK_expected,n_k,n_KK,&
	pf_a(:,z_timeseries(i_sim),:,e_timeseries(i_sim,i_hh)),k_timeseries(2,i_hh))
KK_timeseries(i_sim) = KK_timeseries(i_sim) + k_timeseries(2,i_hh) ! calc total KK
end do
KK_timeseries(i_sim) = KK_timeseries(i_sim)/N_shocks ! unit measure of agents
!write(*,*) "KK_timeseries", KK_timeseries
do i_hh = 1,N_shocks
k_timeseries(1,i_hh) = k_timeseries(2,i_hh)
end do
end do
! step 3: re-estimate parameters
write(*,*) "step"

do i_sim = T_burn+1,T_timeseries
YY(i_sim,1) = log(KK_timeseries(i_sim))
if ( z_timeseries(i_sim) < 1.5 ) then
XX(i_sim,1) = 1
XX(i_sim,2) = log(KK_timeseries(i_sim-1))
XX(i_sim,3) = 0
XX(i_sim,4) = 0
else
XX(i_sim,3) = 1
XX(i_sim,4) = log(KK_timeseries(i_sim-1))
XX(i_sim,1) = 0
XX(i_sim,2) = 0
end if

end do
! just ols stuff here
!iXX = matinv4(matmul(transpose(XX),XX))
call matinv4(matmul(transpose(XX),XX),iXX)
OLS = matmul(matmul(iXX,transpose(XX)),YY)
resid = YY - matmul(XX,OLS)

SSE = 0!transpose(good_resid)*good_resid + transpose(bad_resid)*bad_resid
do i_sim = 1,T_net
SSE = SSE + resid(i_sim,1)**(2)
end do

MSE = 0
meanlogkk = 0
do i_sim= T_burn+1 ,T_timeseries
meanlogkk = meanlogkk + log(KK_timeseries(i_sim))
end do
meanlogkk = meanlogkk/T_net
do i_sim= T_burn+1 ,T_timeseries
MSE = MSE + (log(KK_timeseries(i_sim)) - meanlogkk)**(2) !squared error from just demeaning
end do
Rsq = 1 - SSE/MSE
write(*,*) 'Rsq: ', Rsq

! check for convergence
!do i_KK = 1,2
!OLS(i_KK,1) = good_OLS(i_kk,1)
!OLS(i_KK,2) = bad_OLS(i_kk,1)
!end do
max_diff = maxval(abs(OLS - OLS_lag))
OLS_lag = OLS
a0 = OLS(1,1)
a1 = OLS(2,1)
b0 = OLS(3,1)
b1 = OLS(4,1)
if (max_diff>res_tol) then
write(*,*) 'max diff: ', max_diff
else
converged_outer = 1
end if
end do

! check for convergence
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
200 format(f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,&
	f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x,f25.15,2x) !,f25.15,2x,f25.15,2x

do i_k = 1,n_k
     write(2,200) 
     !, pf_i_apr(i_a,1), pf_i_apr(i_a,2)!grid_a(i_a), pf_c(i_a), pf_a(i_a), pf_v(i_a), pmf(i_a,1), pmf(i_a,2), pmf_init_flat(i_a,1), pmf_init_flat(i_a,2)!, pf_i_apr(i_a,1), pf_i_apr(i_a,2)
end do

return

end subroutine coda


!interpolation code from Phil
subroutine lin_interp_1d(x1,x1i,nx1,pf,o1)

!use params_grid
!use Numerical_Libraries

implicit none

integer, intent(in)             :: nx1
double precision, intent(in)    :: x1(nx1)
double precision, intent(in)    :: x1i
double precision, intent(in)    :: pf(nx1)
double precision, intent(out)   :: o1

double precision :: s1
double precision :: x1i_min
!double precision :: loc1
integer :: loc1

double precision :: xi
double precision :: xi_left
double precision :: xi_right

double precision :: w_2
double precision :: w_1
double precision :: w1(2)

integer :: m1

o1 = 0

s1      = x1(2) - x1(1)
x1i_min = x1i - x1(1)
loc1    = min(nx1 - 1, max(1,floor(x1i_min/s1)+1))

xi = x1i
xi_left = x1(loc1)
xi_right = x1(loc1 + 1)

w_2 = (xi - xi_left)/(xi_right - xi_left)
w_1 = 1 - w_2

w1(1) = w_1
w1(2) = w_2

do m1 = 0,1
    o1 = o1 + w1(m1+1)*pf(loc1+m1)
end do

return

end subroutine lin_interp_1d

subroutine lin_interp_2d(x1,x2,x1i,x2i,nx1,nx2,pf,o1)

!use param_ss_grid
!use Numerical_Libraries

implicit none

integer, intent(in)             :: nx1
integer, intent(in)             :: nx2
double precision, intent(in)    :: x1(nx1)
double precision, intent(in)    :: x2(nx2)
double precision, intent(in)    :: x1i
double precision, intent(in)    :: x2i
double precision, intent(in)    :: pf(nx1,nx2)
double precision, intent(out)   :: o1

double precision :: s1
double precision :: s2
double precision :: x1i_min
double precision :: x2i_min
integer :: loc1
integer :: loc2

double precision :: xi(2)
double precision :: xi_left(2)
double precision :: xi_right(2)

double precision :: w_2(2)
double precision :: w_1(2)
double precision :: w1(2)
double precision :: w2(2)

integer :: m1
integer :: m2

o1 = 0

s1      = x1(2) - x1(1)
x1i_min = x1i - x1(1)
loc1    = min(nx1 - 1, max(1,floor(x1i_min/s1)+1))

s2      = x2(2) - x2(1)
x2i_min = x2i - x2(1)
loc2    = min(nx2 - 1, max(1,floor(x2i_min/s2)+1))

xi(1) = x1i
xi(2) = x2i
xi_left(1) = x1(loc1)
xi_left(2) = x2(loc2)
xi_right(1) = x1(loc1 + 1)
xi_right(2) = x2(loc2 + 1)

w_2(1) = (xi(1) - xi_left(1))/(xi_right(1) - xi_left(1))
w_2(2) = (xi(2) - xi_left(2))/(xi_right(2) - xi_left(2))
w_1(1) = 1 - w_2(1)
w_1(2) = 1 - w_2(2)

w1(1) = w_1(1)
w1(2) = w_2(1)

w2(1) = w_1(2)
w2(2) = w_2(2)

do m1 = 0,1
	do m2 = 0,1
	    o1 = o1 + w1(m1+1)*w2(m2+1)*pf(loc1+m1,loc2+m2)
	end do
end do

return

end subroutine lin_interp_2d


!function matinv4(A) result(B)
subroutine matinv4(A,B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
    double precision, intent(in) :: A(4,4)   !! Matrix
    double precision, intent(out)   :: B(4,4)   !! Inverse matrix
    double precision        :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = &
      1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
       - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
       + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
       - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

    ! Calculate the inverse of the matrix
    B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
    B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
    B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
    B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
    B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
    B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
    B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
    B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
    B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
    B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
    B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
    B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
    B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
    B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
  end subroutine matinv4