all: sim launcher

sim: LAr_dm_B.c LAr_fEk_sec_el.c ran_gen.c LAr_crsec_3.c
	gcc -O3 LAr_dm_B.c LAr_fEk_sec_el.c ran_gen.c LAr_crsec_3.c -lm -openmp -O3 -o sim

launcher: mainLauncher.c
	gcc mainLauncher.c -o launcher

.PHONY: clean outClean

outClean:
	rm -f *.out
	rm -f *.det 

clean:
	rm -f sim
	rm -f launcher