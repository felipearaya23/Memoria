SHELL = /bin/sh
CC = $(CXX)
CFLAGS = -Wall -Wextra -pedantic -Wshadow -funroll-loops -std=c++0x -O3 -DNDEBUG -pthread -march=native
AUX_PAR_FLAGS = -fopenmp
#CFLAGS = -Wall -Wextra -pedantic -Wshadow -std=c++0x -g2 -pthread
#AUX_DISK_FLAGS = -DMONITOR_DISK_USAGE

all: construct_lcp_sequential construct_lcp_parallel

construct_lcp_sequential:
	$(CC) $(CFLAGS) -o construct_lcp_sequential src/em_sparse_phi_src/utils.cpp src/main.cpp $(AUX_DISK_FLAGS)

construct_lcp_parallel:
	$(CC) $(CFLAGS) $(AUX_PAR_FLAGS) -o construct_lcp_parallel src/em_sparse_phi_src/utils.cpp src/main.cpp $(AUX_DISK_FLAGS)

clean:
	/bin/rm -f *.o

nuclear:
	/bin/rm -f construct_lcp_sequential construct_lcp_parallel *.o
