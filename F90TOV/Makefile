F90=gfortran
F90FLAGS=-g -O3

all: mytov mytov_profile TOV_tabular.a

mytov: TOV_tabular.a mytov.F90
	$(F90) $(F90FLAGS) -o mytov mytov.F90 TOV_tabular.a 


mytov_profile: TOV_tabular.a mytov_profile.F90
	$(F90) $(F90FLAGS) -o mytov_profile mytov_profile.F90 TOV_tabular.a 


TOV_tabular.a: TOV_tabular.o 
	ar r TOV_tabular.a TOV_tabular.o 

TOV_tabular.o: TOV_tabular.F90 
	$(F90) $(F90FLAGS) -c TOV_tabular.F90 

clean:
	rm -f *.o *.mod *.a
