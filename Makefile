

all: sim launcher

sim: LAr_dm_B.c LAr_fEk_sec_el.c ran_gen.c LAr_crsec_3.c
	gcc LAr_dm_B.c LAr_fEk_sec_el.c ran_gen.c LAr_crsec_3.c -o sim -lm

launcher: mainLauncher.c
	gcc mainLauncher.c -o launcher

.PHONY: clean

clean:
	rm -f sim
	rm -f launcher