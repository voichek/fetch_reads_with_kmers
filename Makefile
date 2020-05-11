all:fetch_reads.cpp
		g++ -g -O3 fetch_reads.cpp -o fetch_reads -lz 
clean:
		rm -f *.o fetch_reads

