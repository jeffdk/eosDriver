module TOV

  implicit none


  integer,parameter :: N = 20000
  real*8 :: dr = 200.0d0 ! in cm
  real*8 :: TOVphi(N) !nothing
  real*8 :: TOVgmass(N) !in grams
  real*8 :: TOVbmass(N) !in grams
  real*8 :: TOVrho(N) !in grams/cm^3
  real*8 :: TOVpress(N) !in dynes/cm^2
  real*8 :: TOVrad(N) !in cm
  real*8 :: TOVtemp(N) !in MeV
  real*8 :: TOVye(N)
  real*8 :: TOVeps(N) !in ergs/g
  real*8,parameter :: G = 6.6742d-8
  real*8,parameter :: pi = 3.14159265d0
  real*8,parameter :: c = 29979245800.0d0
  real*8,parameter :: precision = 1.0d-12
  real*8 :: eos_rhomax = 0.0d0
  real*8 :: eos_rhomin = 1.0d4
  real*8 :: eos_lrhomin 
  real*8 :: eos_drhoi 
  real*8 :: energy_shift

  real*8,parameter :: rho_gf = 1.61930347d-18
  real*8,parameter :: press_gf = 1.80171810d-39
  real*8,parameter :: eps_gf = 1.11265006d-21
  

  integer :: neos
  integer,parameter :: neos_max = 10000
  real*8  :: eos_lrho(neos_max)
  real*8  :: eos_lprs(neos_max)
  real*8  :: eos_leps(neos_max)

  integer :: last_valid_index
  real*8 :: dummy(15)

contains

  subroutine get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
    
    implicit none

    real*8 :: central_density
    real*8 :: output_rho(N)
    real*8 :: output_mgrav(N)
    real*8 :: output_press(N)

    integer :: last_index
    
  
    !initialize
    call init_TOV

    !do TOV integration
    call do_TOV(central_density)
    
    !set up outputs
    output_rho = TOVrho
    output_mgrav = TOVgmass
    output_press = TOVpress
    last_index = last_valid_index
    
  end subroutine get_TOV_solutions

  subroutine do_TOV(central_density)
      
    implicit none
    
    integer :: i
    integer :: flag_return

    real*8 :: h
    real*8 :: central_density
    real*8 :: temp_rho,temp_temp,temp_ye,temp_eps


    !RK variables
    real*8 :: k1p,k2p,k3p,k4p
    real*8 :: k1m,k2m,k3m,k4m
    real*8 :: k1ph,k2ph,k3ph,k4ph
    real*8 :: k1mb,k2mb,k3mb,k4mb

    flag_return = 0

    !do inner step first
    TOVrad(1) = 0.0d0
    TOVrho(1) = central_density
    TOVgmass(1) = 0.0d0
    TOVbmass(1) = 0.0d0
    TOVphi(1) = 0.0d0

    call tabeos_press_eps(TOVrho(1),TOVpress(1),TOVeps(1))

    do i = 2,N
       TOVrad(i) = dr*dble(i-1)
       h = 1.0d0+TOVeps(i-1)/c**2+TOVpress(i-1)/TOVrho(i-1)/c**2
       call dpdr(TOVgmass(i-1),TOVpress(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1p)
       call dmdr(TOVpress(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1m)
       call dphidr(TOVgmass(i-1),TOVpress(i-1),TOVrho(i-1)*h,TOVrad(i-1),k1ph)
       call dmbdr(TOVgmass(i-1),TOVrho(i-1),TOVrad(i-1),k1mb)

       if (TOVpress(i-1)+dr*0.5d0*k1p.lt.0.0d0) then
          last_valid_index = i-1
!          write(*,*) "here1"
          goto 10
       endif
       
       temp_rho = TOVrho(i-1)
       call get_rho_from_pressure(TOVpress(i-1)+dr*0.5d0*k1p,temp_rho, &
            temp_temp,temp_ye,temp_eps,flag_return)
       if (flag_return.eq.1) then
          last_valid_index = i-1
!          write(*,*) "here2"
          goto 10
       endif
       
       if (temp_rho.lt.2.0d5) then
          last_valid_index = i-1
!          write(*,*) "here3"
          goto 10
       endif

       h = 1.0d0+temp_eps/c**2+(TOVpress(i-1)+0.5d0*dr*k1p)/temp_rho/c**2
       call dpdr(TOVgmass(i-1)+0.5d0*dr*k1m,TOVpress(i-1)+0.5d0*dr*k1p, &
            temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2p)
       call dmdr(TOVpress(i-1)+0.5d0*dr*k1p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2m)
       call dphidr(TOVgmass(i-1)+0.5d0*dr*k1m,TOVpress(i-1)+0.5d0*dr*k1p, &
            temp_rho*h,TOVrad(i-1)+0.5d0*dr,k2ph)
       call dmbdr(TOVgmass(i-1)+0.5d0*dr*k1m,temp_rho,TOVrad(i-1)+0.5d0*dr,k2mb)

       if (TOVpress(i-1)+dr*0.5d0*k2p.lt.0.0d0) then
          last_valid_index = i-1
!          write(*,*) "here4"
          goto 10
       endif
       temp_rho = TOVrho(i-1)
       call get_rho_from_pressure(TOVpress(i-1)+dr*0.5d0*k2p,temp_rho, &
            temp_temp,temp_ye,temp_eps,flag_return)
       if (flag_return.eq.1) then
          last_valid_index = i-1
!          write(*,*) "here5"
          goto 10
       endif
       
       if (temp_rho.lt.2.0d5) then
          last_valid_index = i-1
!          write(*,*) "here6"
          goto 10
       endif

       h = 1.0d0+temp_eps/c**2+(TOVpress(i-1)+dr*0.5d0*k2p)/temp_rho/c**2
       call dpdr(TOVgmass(i-1)+dr*0.5d0*k2m,TOVpress(i-1)+dr*0.5d0*k2p, &
            temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3p)
       call dmdr(TOVpress(i-1)+dr*0.5d0*k2p,temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3m)
       call dphidr(TOVgmass(i-1)+dr*0.5d0*k2m,TOVpress(i-1)+dr*0.5d0*k2p, &
            temp_rho*h,TOVrad(i-1)+0.5d0*dr,k3ph)
       call dmbdr(TOVgmass(i-1)+0.5d0*dr*k2m,temp_rho,TOVrad(i-1)+0.5d0*dr,k3mb)

       if (TOVpress(i-1)+dr*k3p.lt.0.0d0) then
          last_valid_index = i-1
!          write(*,*) "here7"
          goto 10
       endif

       temp_rho = TOVrho(i-1)
       call get_rho_from_pressure(TOVpress(i-1)+dr*k3p,temp_rho, &
            temp_temp,temp_ye,temp_eps,flag_return)
       if (flag_return.eq.1) then
          last_valid_index = i-1
!          write(*,*) "here8"
          goto 10
       endif
       
       if (temp_rho.lt.2.0d5) then
          last_valid_index = i-1
!          write(*,*) "here9"
          goto 10
       endif

       h = 1.0d0+temp_eps/c**2+(TOVpress(i-1)+dr*k3p)/temp_rho/c**2
       call dpdr(TOVgmass(i-1)+dr*k3m,TOVpress(i-1)+dr*k3p, &
            temp_rho*h,TOVrad(i-1)+dr,k4p)
       call dmdr(TOVpress(i-1)+dr*k3p,temp_rho*h,TOVrad(i-1)+dr,k4m)
       call dphidr(TOVgmass(i-1)+dr*k3m,TOVpress(i-1)+dr*k3p, &
            temp_rho*h,TOVrad(i-1)+dr,k4ph)
       call dmbdr(TOVgmass(i-1)+dr*k3m,temp_rho,TOVrad(i-1)+dr,k4mb)
       
       TOVpress(i) = TOVpress(i-1) + dr/6.0d0*(k1p + 2.0d0*k2p + &
            2.0d0*k3p + k4p)
       if (TOVpress(i).lt.0.0d0) then
          last_valid_index = i-1
!          write(*,*) "here10"
          goto 10
       endif
       TOVgmass(i) = TOVgmass(i-1) + dr/6.0d0*(k1m + 2.0d0*k2m + &
            2.0d0*k3m + k4m)
       TOVphi(i) = TOVphi(i-1) + dr/6.0d0*(k1ph + 2.0d0*k2ph + &
            2.0d0*k3ph + k4ph)
       TOVbmass(i) = TOVbmass(i-1) + dr/6.0d0*(k1mb + 2.0d0*k2mb + &
            2.0d0*k3mb + k4mb)
       
       temp_rho = TOVrho(i-1)
       call get_rho_from_pressure(TOVpress(i),temp_rho, &
            temp_temp,temp_ye,temp_eps,flag_return)
       if (flag_return.eq.1) then
          last_valid_index = i-1
!          write(*,*) "here11"
          goto 10
       endif

       TOVrho(i) = temp_rho

       call tabeos_press_eps(TOVrho(i),TOVpress(i),TOVeps(i))


       if (TOVpress(i)/TOVpress(1).lt.1.0d-10) then
          last_valid_index = i
!          write(*,*) "here12"
          goto 10
       endif

    enddo
    
    
    stop "Consider making the TOV radial space bigger, you are running out of space for a full solution"

10  continue

  end subroutine do_TOV
  
  subroutine get_rho_from_pressure(xpress,xrho,temp_temp,temp_ye,xeps,flag_return)
    implicit none
    
    real*8 :: xpress
    real*8 :: xrho,xeps,leps
    real*8 :: temp_temp,temp_ye,temp_eps
    real*8 :: rho_guess, rho_guess2
    real*8 :: press_guess,press_guess2
    real*8 :: mydpdrho
    real*8 :: fac 
    
    integer :: counter
    integer :: flag_return

    logical :: cont

    real*8  :: lprec

    lprec = precision

    fac = 1.0d0
    rho_guess = xrho
    cont = .true.
    counter = 0
!    write(6,*) "****************"
    do while(cont)
       counter = counter + 1
       rho_guess2 = rho_guess*1.0001d0
!       write(*,*) counter, rho_guess, xpress,press_guess

       call tabeos_press_eps(rho_guess2,press_guess2,leps)
       call tabeos_press_eps(rho_guess,press_guess,leps)

       mydpdrho = (press_guess2-press_guess)/(rho_guess2-rho_guess)
       if (mydpdrho.lt.0.0d0) then
!          write(6,"(A30,1P1E15.6,i10)") "table is crap, dpdrho.lt.0",rho_guess,counter
       endif
       
       if (abs(1.0d0-press_guess/xpress).lt.lprec) then
          cont = .false.
          xrho = rho_guess
          xeps = leps
       else
          if (fac*(xpress-press_guess)/mydpdrho/rho_guess.gt.0.1d0) then
             rho_guess = 0.99d0*rho_guess
          else
             rho_guess = rho_guess + fac*(xpress-press_guess)/mydpdrho
          endif
       endif
       
       if (counter.gt.100) then
          fac = 0.01d0
       endif

       if (counter.gt.1000) then
          fac = 0.001d0
          lprec = lprec*100
       endif

       if (counter.gt.10000) then
          lprec = lprec*1000
       endif

!       write(6,"(i8,1P10E15.6)") counter,abs(1.0d0-press_guess/xpress),&
!            press_guess,xpress,rho_guess, (xpress-press_guess)/mydpdrho, precision
       
       if (rho_guess.le.eos_rhomin) then
          cont = .false.
          xrho = eos_rhomin
          flag_return = 1
       endif

       if(counter.gt.50000) then
          write(6,*) "problem in rho(press)"
          write(6,"(1P10E15.6)") rho_guess,abs(1.0d0-press_guess/xpress),lprec,eos_rhomin
          stop "fucked up"
       endif
       
    enddo
    
  end subroutine get_rho_from_pressure

  subroutine init_TOV
    
    implicit none

!!$    allocate(TOVphi(N))
!!$    allocate(TOVgmass(N))    
!!$    allocate(TOVbmass(N))    
!!$    allocate(TOVrho(N))
!!$    allocate(TOVpress(N))
!!$    allocate(TOVrad(N))
!!$    allocate(TOVtemp(N))
!!$    allocate(TOVye(N))
!!$    allocate(TOVeps(N))

    TOVphi(:) = 0.0d0
    TOVgmass(:) = 0.0d0
    TOVbmass(:) = 0.0d0
    TOVrho(:) = 0.0d0
    TOVpress(:) = 0.0d0
    TOVrad(:) = 0.0d0
    TOVtemp(:) = 0.0d0
    TOVye(:) = 0.0d0
    TOVeps(:) = 0.0d0
    

  end subroutine init_TOV

!!$  subroutine deallocate_TOV
!!$
!!$    implicit none
!!$
!!$    deallocate(TOVphi)
!!$    deallocate(TOVgmass)
!!$    deallocate(TOVbmass)
!!$    deallocate(TOVrho)
!!$    deallocate(TOVpress)
!!$    deallocate(TOVrad)
!!$    deallocate(TOVtemp)
!!$    deallocate(TOVye)
!!$    deallocate(TOVeps)
!!$
!!$  end subroutine deallocate_TOV
  
  subroutine dpdr(mass,press,rho,rad,xx)
    
    implicit none

    !inputs
    real*8 mass,press,rho,rad
    
    !outputs
    real*8 xx

    if (mass.eq.0.0d0) then
       xx = -G*rho*4.0d0*pi*press*rad/c**2
    else
       xx = -G*rho*mass/rad**2*(1.0d0+4.0d0*pi*press/c/c*rad**3/mass)* &
            (1.0d0-2.0d0*G*mass/c/c/rad)**(-1)
    endif

  end subroutine dpdr

  subroutine dmdr(press,rho,rad,xx)
    
    implicit none

    !inputs
    real*8 mass,press,rho,rad

    !Output
    real*8 xx

    xx = 4.0d0*pi*rad**2*(rho-press/c**2)
    
  end subroutine dmdr
  
  subroutine dmbdr(mass,rho,rad,xx)
    
    implicit none
    
    !inputs
    real*8 mass,press,rho,rad
    
    !Output
    real*8 xx
    
    if (mass.eq.0.0d0) then
       xx = 4.0d0*pi*rad**2*rho
    else
       xx = 4.0d0*pi*rad**2*rho*(1.0d0-2.0d0*G*mass/(rad*c**2))**(-0.5d0)
    endif
    
  end subroutine dmbdr

  subroutine dphidr(mass,press,rho,rad,xx)

    implicit none
    
    !inputs
    real*8 mass,press,rho,rad
    
    !outputs
    real*8 xx
    
    if (mass.eq.0.0d0) then
       xx = G*4.0d0*pi*press*rad/c**4
    else
       xx = (G/c**2)*mass/rad**2*(1.0d0+4.0d0*pi*press/c/c*rad**3/mass)* &
            (1.0d0-2.0d0*G*mass/c/c/rad)**(-1)
    endif
    
  end subroutine dphidr

  subroutine readtable(tablename)

    implicit none

    real*8  :: buffer1,buffer2,buffer3
    integer :: iostatus
    integer :: count
    character(*) :: tablename

    count = 0
    open(unit=666,file=trim(adjustl(tablename)),action="read")
    read(666,*) energy_shift
    do
       read(666,*,iostat=iostatus) buffer1,buffer2,buffer3
       if(iostatus .ne. 0) exit
       count = count + 1
       if(count > neos_max) then
          stop "Too many entries in EOS table!"
       endif
       eos_lrho(count) = dlog10(10.0d0**buffer1)
       eos_lprs(count) = dlog10(10.0d0**buffer2 / press_gf)
       eos_leps(count) = dlog10(10.0d0**buffer3 / eps_gf)
    enddo
    close(666)
    neos = count
    eos_rhomax = 10.0d0**(eos_lrho(neos))
    eos_rhomin = 10.0d0**(eos_lrho(1))
    eos_lrhomin = dlog10(eos_rhomin)
    eos_drhoi  = 1.0d0/(eos_lrho(2) - eos_lrho(1))
    
    write(6,*) "Read EOS table ", trim(adjustl(tablename))
    write(6,"(A10,i7,A10,1P10E15.6)") "nrho: ",neos," rho_min/max: ",eos_rhomin,eos_rhomax
    write(6,"(A10,1P10E15.6)") "eps shift", energy_shift

  end subroutine readtable

  subroutine tabeos_press_eps(rho,press,eps)

    implicit none
    real*8,intent(IN) :: rho
    real*8,intent(OUT) :: press, eps
    real*8 :: lrho, lprs, leps
    integer :: irho

    lrho = dlog10(rho)
    irho = 2 + int( (lrho - eos_lrhomin - 1.0d-10) * eos_drhoi )

    if(irho < 2) then
       write(6,"(i8,1P10E15.6)") irho,lrho,rho
       stop "outside table bounds!"
    endif
    if(irho.gt.neos) stop "rho > rho_max"

    lprs = (eos_lprs(irho) - eos_lprs(irho-1)) * eos_drhoi * &
         (lrho - eos_lrho(irho-1)) + eos_lprs(irho-1)

    leps = (eos_leps(irho) - eos_leps(irho-1)) * eos_drhoi * &
         (lrho - eos_lrho(irho-1)) + eos_leps(irho-1)

    press = 10.0d0**lprs
    eps = 10.0d0**leps - energy_shift
    
  end subroutine tabeos_press_eps



end module TOV
