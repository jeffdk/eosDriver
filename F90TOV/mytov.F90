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

  integer, parameter :: nrhos = 400
  real*8 :: rhos(nrhos),mgravs(nrhos),mbarys(nrhos),rads(nrhos)
  real*8 :: rho_min = 6.0d14
  real*8 :: rho_max = 7.5d15
  real*8 :: dlrho, lrho
  real*8 :: buff1,buff2

  character(len=512) :: outfilename
  character(len=512) :: eosfilename

  allocate(output_rho(N))
  allocate(output_mgrav(N))
  allocate(output_press(N))

  eosfilename = "eostable_ye=00.150_s=00.500.dat"
  outfilename = "tov_sequence_"//trim(adjustl(eosfilename))

  write(6,*) outfilename

  call readtable(trim(adjustl(eosfilename)))

  dlrho = (log10(rho_max) - log10(rho_min)) / (nrhos - 1)
  do i=1,nrhos
     rhos(i) = 10.0d0**(log10(rho_min) + (i-1)*dlrho)
  enddo

  do i=1,nrhos
     central_density = rhos(i)
     call get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
     mgravs(i) = output_mgrav(last_index)/msun
     mbarys(i) = TOVbmass(last_index)/msun
     rads(i) = TOVrad(last_index)
     write(6,"(i5,1P1E15.6,1P3G12.4)") i,central_density,mgravs(i),mbarys(i),rads(i)/1.0e5
  enddo

  open(666,file=outfilename,action="write")
  do i=1,nrhos
     write(666,"(i5,1P10E15.6)") i,rhos(i),mgravs(i),mbarys(i),rads(i)
  enddo
  close(666)


end program mytov

#if 0
  central_density = 1.0d15
  call get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
  write(6,*) output_mgrav(last_index)/msun, last_index

  open(unit=666,file="outfile.dat")
  do i=1,last_index
     write(666,"(1P10E18.9)") TOVrad(i),TOVrho(i),TOVpress(i),TOVgmass(i)
  enddo
  close(666)

  write(6,*) "******************************************"
  write(6,*) "******************************************"
  write(6,*) "******************************************"

  central_density = 1.00001d15
  call get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
  write(6,*) output_mgrav(last_index)/msun, last_index

  open(unit=666,file="outfile2.dat")
  do i=1,last_index
     write(666,"(1P10E18.9)") TOVrad(i),TOVrho(i),TOVpress(i),TOVgmass(i)
  enddo

  central_density = 1.0d15
  call get_TOV_solutions(central_density,output_rho,output_press,output_mgrav,last_index)
  write(6,*) output_mgrav(last_index)/msun, last_index

#endif
