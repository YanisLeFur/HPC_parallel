CXX=nvcc
LD=${CXX}
CXXFLAGS+=-Xcompiler="-Wall -Wextra -Werror" -std c++11 -O3
LDFLAGS+=-lm $(CXXFLAGS)

OBJS_NBODY=nbody.o reader.o

all: nbody-code-entry nbody-code-row


nbody-code-row: $(OBJS_NBODY) nbody_bruteforce_row.o
	$(LD) -o $@ $(OBJS_NBODY) nbody_bruteforce_row.o $(LDFLAGS) 


nbody-code-entry : $(OBJS_NBODY) nbody_bruteforce_entry.o
	$(LD) -o $@ $(OBJS_NBODY) nbody_bruteforce_entry.o $(LDFLAGS) 

nbody_bruteforce_row.o: nbody_bruteforce.cu
	$(CXX) $(CXXFLAGS) -DPER_ROW -c $< -o $@

nbody_bruteforce_entry.o: nbody_bruteforce.cu
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -Rf nbody-code *.o *~  gmon.out


# CXX=nvcc
# LD=${CXX}
# CXXFLAGS+=-Xcompiler="-Wall -Wextra -Werror" -std c++11 -O3
# LDFLAGS+=-lm $(CXXFLAGS)

# OBJS_NBODY=nbody.o reader.o

# all: clean nbody-code

# nbody-code: $(OBJS_NBODY) nbody_bruteforce.o
# 	$(LD) -o $@ $(OBJS_NBODY) nbody_bruteforce.o $(LDFLAGS) 

# nbody_bruteforce.o: nbody_bruteforce.cu
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean:
# 	rm -Rf nbody-code *.o *~  gmon.out
