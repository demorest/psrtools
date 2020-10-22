PROGS = autotoa normalize_rms
all: $(PROGS)
install: $(PROGS)
	cp $(PROGS) $(PREFIX)/bin

autotoa: autotoa.C
	$(CXX) $(CXXFLAGS) `psrchive --cflags` -o autotoa autotoa.C $(LDFLAGS) `psrchive --libs` 

normalize_rms: normalize_rms.C
	$(CXX) $(CXXFLAGS) `psrchive --cflags` -o normalize_rms normalize_rms.C $(LDFLAGS) `psrchive --libs` 

clean:
	rm -f $(PROGS) *.o
