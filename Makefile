PROGS = autotoa normalize_rms
all: $(PROGS)
install: $(PROGS)
	cp $(PROGS) $(PSR)/bin

autotoa: autotoa.C
	$(CXX) -g `psrchive --cflags` -o autotoa autotoa.C `psrchive --libs` 

parallactic: parallactic.C
	$(CXX) -g `psrchive --cflags` -o parallactic parallactic.C `psrchive --libs` 

check_model_phase: check_model_phase.C
	$(CXX) -g `psrchive --cflags` -o check_model_phase check_model_phase.C `psrchive --libs` 

normalize_rms: normalize_rms.C
	$(CXX) -g `psrchive --cflags` -o normalize_rms normalize_rms.C `psrchive --libs` 

dm_stat_search: dm_stat_search.C
	$(CXX) -O2 `dspsr_cflags` `psrchive --cflags` -o $@ $< `dspsr_ldflags`

clean:
	rm -f $(PROGS) *.o
