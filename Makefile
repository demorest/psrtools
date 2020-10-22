PROGS = autotoa normalize_rms
all: $(PROGS)
install: $(PROGS)
	cp $(PROGS) $(PREFIX)/bin

autotoa: autotoa.C
	$(CXX) -g `psrchive --cflags` -o autotoa autotoa.C `psrchive --libs` 

normalize_rms: normalize_rms.C
	$(CXX) -g `psrchive --cflags` -o normalize_rms normalize_rms.C `psrchive --libs` 

clean:
	rm -f $(PROGS) *.o
