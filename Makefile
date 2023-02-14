# TODO: Change Makefile back to original!
# TODO: Make sure to use -O3 when answering any of the report questions.

ifdef USE_INT
MACRO = -DUSE_INT
endif

#compiler setup
CXX = g++
WFLAGS = -Wall -Werror -Wno-error=unknown-pragmas
CXXFLAGS = -std=c++14 -g -pthread $(WFLAGS) $(MACRO)

COMMON= core/utils.h core/cxxopts.h core/get_time.h core/graph.h core/quick_sort.h
SERIAL= triangle_counting page_rank
PARALLEL= triangle_counting_parallel page_rank_parallel
ALL= $(SERIAL) # $(PARALLEL)

all : $(ALL)

% : %.cpp $(COMMON)
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)