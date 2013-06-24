program mytov
  use tov
  implicit none
  
  real*8 :: central_density 
  real*8,parameter :: msun = 1.9891d33
  real*8,allocatable :: output_rho(:)
  real*8,allocatable :: output_mgrav(:)
  real*8,allocatable :: output_press(:)
  integer :: last_index
  integer :: i

  integer, parameter :: nrhos = 1000
  real*8 :: rhos(nrhos),mgravs(nrhos),mbarys(nrhos),rads(nrhos)
  real*8 :: rho_min = 7.0d14
  real*8 :: rho_max = 9.0d14
  real*8 :: dlrho, lrho
  real*8 :: buff1,buff2

  character(len=512) :: outfilename
  character(len=128) :: eosfilename


  allocate(output_rho(N))
  allocate(output_mgrav(N))
  allocate(output_press(N))

  call getarg(1, eosfilename)
!  eosfilename = "eostable_ye=00.150_s=00.500.dat"
  outfilename = "tov_profile_135Msun_"//trim(adjustl(eosfilename))

  write(6,*) eosfilename
  write(6,*) outfilename

  call readtable(trim(adjustl(eosfilename)))

  central_density = 5.578e14
  call get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
  write(6,*) output_mgrav(last_index)/msun, last_index

  open(unit=666,file=trim(adjustl(outfilename)))
  do i=1,last_index
     write(666,"(1P10E18.9)") TOVrad(i),TOVrho(i),TOVpress(i),TOVgmass(i)
  enddo
  close(666)


end program mytov

