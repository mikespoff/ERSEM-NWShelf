!-----------------------------------------------------------------------
! Copyright 2014 Plymouth Marine Laboratory
!
! This file is part of the SSB-ERSEM library.
!
! SSB-ERSEM is free software: you can redistribute it and/or modify it 
! under the terms of the Lesser GNU General Public License as published 
! by the Free Software Foundation, either version 3 of the License, or 
! (at your option) any later version.
!
! SSB-ERSEM is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser 
! GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License 
! along with SSB-ERSEM. If not, see <http://www.gnu.org/licenses/>.
!
! Address:
! Plymouth Marine Laboratory
! Prospect Place, The Hoe
! Plymouth, PL1 3DH, UK
!
! Email:
! ssm@pml.ac.uk
!-----------------------------------------------------------------------
#include"ppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: bio_ersem
!
! !DESCRIPTION:
!  TODO - Descripton
!\\
!\\
! !INTERFACE:   
   module bio_ersem
!
! Author: M. Butenschön, PML - August 2008
!
! !USES:
   use ersem_constants, only: fp8, seconds_per_day
   use ersem_variables, only: eLogU, eDbgU
   use pelagic_variables, only: n_comp, i_state, ccc, sccc, itrXccc, EIR
   use benthic_variables, only: n_compben, i_stateben, ccb, sccb, itrXccb
   use ersem_variables, only: ncdfInstOut, ncdfDailyOut, &
                              ncdfWeeklyOut, ncdfMonthlyOut
   use meanflow, only: bioshade,h,w
!
   implicit none
!
!  Default all is private
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public step_ersem,gotm2ersem,getSecsOfMonth,getWeekday
!
! !PUBLIC DATA MEMBERS:
   logical, public   :: bio_calc=.true.,splitting=.true. &
           ,bioshade_feedback=.true.
   integer, public    :: ncc,split_factor=1,nbudget, &
                         w_adv_discr=6,w_adv_ctr=0,ode_method=1 
   real(fp8), public, dimension(:,:), allocatable :: cc,ccben,ccpel,sccpel
   real(fp8), public     :: cnpar=0.5,kc=0.03
   integer, public :: nmon,ny,nsmonth,nmonth0, time_second
!
! !PRIVATE DATA MEMBERS:
!  logical                   :: bio_eulerian=.true.
!  integer                   :: bio_npar=10000
!
! !REVISION HISTORY:
!  Original author(s):
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: step_ersem
!
! !DESCRIPTION:
!  TODO - description
!\\
!\\
! !INTERFACE:
   subroutine step_ersem(nlev,dt,nnow)
!
! !USES:
!  From ersem:
   use budget, only: calc_budget,nmonthly,nwoff,ndoff,nmoff
   use ncdfersem, only: write_to_ncdfersem,ncdfI,ncdfD,ncdfW,ncdfM
!  From gotm:
   use util, only:Dirichlet,Neumann,flux
   use turbulence, only:nuh
   use output, only:ts
   use airsea, only: I_0
   use time, only: MaxN
#ifdef PPNCDFRESTART
   use ncdfRestartErsem, only:ncdfR,ncdfWriteRestart, &
        ncdfCloseRestart,ncdfR,ccrst,cbrst,writeErsemRestart
#endif
   use reset, only: reset_ersem
   use budget, only: set_fixed_quota
!
! !INPUT PARAMETERS:
   integer,  intent(in)                :: nlev,nnow
   real(fp8), intent(in)                :: dt
!
! !LOCAL VARIABLES:
   real(fp8)                  :: RelaxTau(0:nlev)
   integer :: j,jj,i,nsec
   real(fp8),dimension(0:nlev) :: Lsour,conc,Qsour
   integer, parameter        :: adv_mode_0=0
   integer, parameter        :: adv_mode_1=1
   real(fp8)                  :: dt_eff
   integer                   :: split
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
! 
   if (bio_calc) then

      RelaxTau = 1.e15
      Lsour=0.e0

      if (splitting) then
         sccpel=0.e0
      else
         ! call ERSEM (sets sccc->sccpel,sccb):
         call do_ersem_rates(nlev)
      end if

!     Advectin-diffusion and time integration:
      do j=1,i_state
         if(itrXccc(j).gt.10) then
            Qsour=sccpel(j,:)
            conc=ccpel(j,:)
!           do advection step due to settling or rising
!            call adv_center(nlev,dt,h,h,ws(j,:),flux,                   &
!                 flux,_ZERO_,_ZERO_,w_adv_discr,adv_mode_1,cc(j,:))

!           do advection step due to vertical velocity
            if(w_adv_ctr .ne. 0) then
               call adv_center(nlev,dt,h,h,w,flux,                   &
                    flux,0._fp8,0._fp8,w_adv_ctr,adv_mode_0,conc)
            end if
            
!           do diffusion step
            call diff_center(nlev,dt,cnpar,1,h,Neumann,Neumann,&
                0._fp8,0._fp8,nuh,Lsour,Qsour,RelaxTau,conc,conc)
            !    0.0d0,0.0d0,nuh,Lsour,sccpel(j,:),RelaxTau,ccpel(j,:),ccpel(j,:))
            ccpel(j,:)=conc
         end if
      end do

      ! Pass updated pelagic states to ERSEM:
      do i=1,I_STATE
         ccc(1:,i)=ccpel(i,nlev:1:-1)
      end do

      if (splitting) then
         ! Pass transported state into full array
         cc(1,1:N_COMP*I_STATE)=reshape(ccpel(:,1:),(/N_COMP*I_STATE/))
         do split=1,split_factor
            dt_eff=dt/float(split_factor)
   
!           Very important for 3D models to save extra 3D field:
            bioshade=1
            ! Call ERSEM (computes rates(sccc,sccb)->pp) and solve for cc
            ! Mind, that the patankar schemes do not work for ERSEM as 
            ! rates are not split in production and destruction part!
            call ode_solver(ode_method,1,ncc,dt_eff,cc,do_ersem_all,do_ersem_ppdd_dummy)
            ! extract pelagic states for advection-diffusion (and ERSEMS ccc):
            ccpel(:,1:)=reshape(cc(1,1:N_COMP*I_STATE),(/I_STATE,N_COMP/))
            ! pass updated pelagic states to ERSEM:
            do i=1,I_STATE
              ccc(1:,i)=ccpel(i,nlev:1:-1)
            end do
            ! pass benthic states to ERSEM
            ccb(1,:)=cc(1,N_COMP*I_STATE+1:)
         end do
      else
         ! solving method fixed to match the one of the advection-diffusion
         ! operator for mass conservation
         ! (The patankar schemes do not work for ERSEM as rates are not 
         ! split  in production and destruction part)
         call ode_solver(1,1,I_STATEBEN,dt,ccben,get_ersem_benthos,get_ersem_benthos_dummy)
         ! pass updated benthic states to ERSEM:
         ccb(1,:)=ccben(1,1:I_STATEBEN)
      end if

      call bioshade_ersem(nlev,bioshade_feedback)

      nsec=nint(nnow*dt)
      ! add seconds to public data member
      time_second=nsec
      if((nsec-nmonth0).eq.nsmonth) then 
          nmoff=nsmonth/2
          nmonthly=0
          nmonth0=nsec
          call getSecsOfMonth(nmon,ny,nsmonth)
      end if

      call set_fixed_quota
      call calc_budget(nnow,dt/86400._fp8,nbudget)

      ! Write ERSEM output:
      if (ncdfInstOut .and. mod(nnow,nbudget).eq.0) then
          write(6,*) 'ERSEM output at ',nsec,'s'
          call write_to_ncdfersem(ncdfI,nsec)
      end if

      if (ncdfDailyOut .and. MOD(nsec+ndoff,86400).eq.0) call write_to_ncdfersem(ncdfD,nsec-43200)
      if (ncdfWeeklyOut .and. MOD(nsec+nwoff,604800).eq.0) call write_to_ncdfersem(ncdfW,nsec-302400)
      if (ncdfMonthlyOut .and. nmonthly.eq.1) call write_to_ncdfersem(ncdfM,nsec-nmoff)

      call reset_ersem

   end if

#ifdef PPNCDFRESTART
   if (writeErsemRestart .and. nnow.eq.MaxN) then
      ccrst=>ccc
      cbrst=>ccb
      call ncdfWriteRestart(ncdfR,ccrst,cbrst)
   endif
#endif

   return

   end subroutine step_ersem
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_ersem_rates
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE
   subroutine do_ersem_rates(nlev)
!
! !USES:
   use ersem, only: ersem_loop
!
! !INPUT PARAMETERS:
   integer,intent(in) :: nlev
!
! !LOCAL VARIABLES:
   integer :: i, j
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
! 
   ! Get environmental forcing:
   call gotm2ersem(nlev)

   ! Ersem main loop
   call ersem_loop

! #MDS Forced addition of nutrients
! phosphate -> sccc(1,3)
! nitrate -> sccc(1,4)
! ammonium -> sccc(1,5)
! silicate -> sccc(1,6)
! 1 day = 86400 secs
! YEARS: 365; 730; 1095; 1461; 1826; 2191; 2556; 2922

! major injection of nutrients
   if (time_second .ge. 86400*(1461+20)) then
      if (time_second .le. 86400*(1461+21)) then
         sccc(1,4) = sccc(1,4) + 500
      endif
   endif

! slow influx of nutrients
   if (time_second .ge. 86400*30) then
  	do i=1,10
   	     sccc(1,3) = sccc(1,3) + 0.002
   	     sccc(1,6) = sccc(1,6) + 0.015
   	end do
   end if

   if (time_second .ge. 86400*800) then
  	do i=1,10
   	     sccc(1,4) = sccc(1,4) + 0.03
   	end do
   end if

  
! End addition of nitrate

   ! Pass ersem rates to gotm:
   do i=1,I_STATE
        sccpel(i,1:)=sccc(N_COMP:1:-1,i)/seconds_per_day 
   end do

   end subroutine do_ersem_rates
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ersem_benthos
!
! !DESCRIPTION:
!  TODO - description
!\\
!\\
! !INTERFACE:
   subroutine get_ersem_benthos(first,nben,nlev,cc,rhs)
!
! !USES:
!
! !INPUT PARAMETERS:
   logical, intent(in) :: first
   integer, intent(in) :: nben,nlev
   real(fp8),intent(in) :: cc(1:nben,0:nlev)
   real(fp8),intent(out) :: rhs(1:nben,0:nlev)
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!  
   rhs=0.0e0
   rhs(1,1:)=sccb(1,1:)/seconds_per_day

   end subroutine get_ersem_benthos
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ersem_benthos_dummy
!
! !DESCRIPTION:
!  TODO - description
!\\
!\\
! !INTERFACE:
   subroutine get_ersem_benthos_dummy(first,nben,nlev,cc,pp,dd)
!
! !USES:
!
! !INPUT PARAMETERS:
   logical, intent(in) :: first
   integer, intent(in) :: nben,nlev
   real(fp8),intent(in) :: cc(1:nben,0:nlev)
   real(fp8),intent(out) :: pp(1:nben,1:nben,0:nlev),dd(1:nben,1:nben,0:nlev)
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
! 
   PPWRITELOG "get_ersem_benthos_dummy ERROR: pp-dd source sink separation is not implemented!!!" 
   PPWRITELOG "This subroutine should not be called."
   PPWRITELOG "Use a non-patankar integration scheme."

   stop

   end subroutine get_ersem_benthos_dummy
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_ersem_ppdd_dummy
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
   subroutine do_ersem_ppdd_dummy(first,ndum,ntot,cc,pp,dd)
!
! !USES:
!
! !INPUT PARAMETERS:
   logical, intent(in) :: first
   integer, intent(in) :: ndum,ntot
   real(fp8),intent(in) :: cc(1:ndum,0:ntot)
   real(fp8),intent(out) :: pp(1:ndum,1:ndum,0:ntot),dd(1:ndum,1:ndum,0:ntot)
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!   
   PPWRITELOG "ERROR: pp-dd source sink separation is not implemented!!!" 
   PPWRITELOG "This subroutine should not be called."
   PPWRITELOG "Use a non-patankar integration scheme."

   stop

   end subroutine do_ersem_ppdd_dummy
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_ersem_all
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
   subroutine do_ersem_all(first,ndum,ntot,cc,rhs)
!
! !USES:
   use ersem, only: ersem_loop
!
! !INPUT PARAMETERS:
   logical, intent(in) :: first
   integer, intent(in) :: ndum,ntot
   real(fp8),intent(in) :: cc(1:ndum,0:ntot)
   real(fp8),intent(out) :: rhs(1:ndum,0:ntot)
!
! !LOCAL VARIABLES:
   integer :: i
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   rhs=0.e0

   ! Get environmental forcing:
   call gotm2ersem(N_COMP)

   ! Ersem loop
   call ersem_loop

   ! Pass ersem rates to gotm:
   do i=1,I_STATE
      sccpel(i,1:)=sccc(N_COMP:1:-1,i)/seconds_per_day 
   end do

   rhs(1,1:N_COMP*I_STATE)=reshape(sccpel(:,1:),(/N_COMP*I_STATE/))
   rhs(1,N_COMP*I_STATE+1:)=sccb(1,:)/seconds_per_day

   end subroutine do_ersem_all
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: bioshade_ersem
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
   subroutine bioshade_ersem(nlev,bioshade_feedback)
!
! !USES:
   use pelagic_variables, only: P1c,P2c,P3c,P4c,R4c,R6c,R8c, &
                                ESS,EPS0X,xEPS
   use light_extinction, only: EPSR6X,EPSESSX
!
! !INPUT PARAMETERS:
   integer, intent(in)   :: nlev
   logical, intent(in)   :: bioshade_feedback
!
! !LOCAL VARIABLES:
   real(fp8) ::  add
   integer               :: i,j
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!   
   add = 0.e0
   j=0
   do i=nlev,1,-1
      j=j+1
      !add=add+h(i)*(P1c(j)+P2c(j)+P3c(j)+P4c(j)+ &
      !          R4c(j)+R6c(j)+R8c(j))
      !if (bioshade_feedback) bioshade(i)=exp(-kc*add)
      add=add*h(i)*(xEPS(j)-EPS0X(j)-EPSR6X*(R6c(j)+R4c(j) &
        +R8c(j))-EPSESSX*ESS(j))
      if (bioshade_feedback) bioshade(i)=exp(-kc*add)
   end do

   end subroutine bioshade_ersem
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: gotm2ersem \label{sec:pelbac}
!
! !DESCRIPTION:
!  TODO - description
!\\
!\\
! !INTERFACE:
   subroutine gotm2ersem(nlev)
!
! !USES:
   use pelagic_variables, only: ESS,WND,ETW,x1x,EIR,SurfaceEIR,xEPS,erw, &
                                pdepth,pvol,parea,BoxDepth,BoxFaces
   use meanflow, only:t,s,z,h,buoy,rho_0,gravity
   use airsea, only:I_0,w
!   use budget, only:nuv
   use turbulence, only:nuh
   use light_extinction, only: calculate_extinction
!
! !INPUT PARAMETERS:
   integer,intent(in) :: nlev
!
! !LOCAL VARIABLES:
   integer :: i
   real(fp8) :: buffer,xtnc
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   BoxDepth(1,1,:)=z(nlev:1:-1)
   BoxFaces(1,1,1:nlev)=z(nlev:1:-1)+.5*h(nlev:1:-1)
   BoxFaces(1,1,nlev+1)=z(1)-.5*h(1)
   WND=w
   ETW=t(nlev:1:-1)
   x1x=s(nlev:1:-1)
   pdepth=h(nlev:1:-1)
   pvol=parea*pdepth
   call calculateESS
   call calculate_extinction
   SurfaceEIR = I_0
   buffer=I_0
   do i=1,N_COMP
      xtnc=xEPS(i)*pdepth(i)
      EIR(i)=buffer/xtnc*(1-exp(-xtnc))
      buffer=buffer*exp(-xtnc)
   end do  
!   nuv=.5*(nuh(nlev:1:-1)+nuh(nlev-1:0:-1))
   erw=-buoy(nlev:1:-1)*rho_0/gravity+rho_0

   return

  end subroutine gotm2ersem
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: calculateESS
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
   subroutine calculateESS
!
! !USES:
!   use sediment, only: C
   use pelagic_variables, only : ess,essXR
!
! !LOCAL VARIABLES:
   integer :: i
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   do i=1,n_comp
      ess(i)=essXR
   end do

   return

   end subroutine calculateESS
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getSecsOfMonth
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
   subroutine getSecsOfMonth(nmonth,nyear,nseconds)
!
! !USES:
   use time, only: time_diff,julian_day
!
! !INPUT PARAMETERS:
   integer :: nmonth,nyear
   integer :: nseconds
!
! !LOCAL VARIABLES:
   integer :: jd1,jd2
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   call julian_day(nyear,nmonth,1,jd1) 
   if (nmonth.eq.12) then
      nmonth=1
      nyear=nyear+1
   else
      nmonth=nmonth+1
   end if
   call julian_day(nyear,nmonth,1,jd2) 
   nseconds=time_diff(jd2,0,jd1,0)

   end subroutine getSecsOfMonth
!
!EOC
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: getWeekday
!
! !DESCRIPTION:
!  Calculates the day of the week from a calendar date
!  Sunday=0 .... Friday=6
!\\
!\\
! !INTERFACE:
   subroutine getWeekday(nyear,nm,nd,nwd)
!
! !USES:
   use time, only: time_diff,julian_day

! !INPUT PARAMETERS:
   integer,intent(in) :: nyear,nm,nd
!
! !OUTPUT PARAMETERS:
   integer,intent(out) :: nwd
!
! !LOCAL VARIABLES:
   integer :: ndd,jdd,jd,rd
!
!EOP
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOC
!
   ! Adoption of Gregorian calendar in UK and US after 02/09/1752
   call julian_day(1752,9,2,rd)
   call julian_day(nyear,nm,nd,jd) 

   ! Get Doomsday of current year:
   if (jd.gt.rd) then
      ndd=MOD(2+nyear+nyear/4-nyear/100+nyear/400,7)
   else
      ndd=MOD(nyear+nyear/4,7)
   end if
   call julian_day(nyear,4,4,jdd) ! reference to Doomsday April 4
   call julian_day(nyear,nm,nd,jd) 

   ! Get current day of the week:
   nwd=MOD(364+ndd+jd-jdd,7)
   write(6,*) 'Weekday:',nwd

   end subroutine getWeekday
!
!EOC
!-----------------------------------------------------------------------

   end module bio_ersem

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
