DEBUG=n

ifeq ($(DEBUG),n)
	OPTFLAG= -O3
else
	OPTFLAG= -O0 -g -DDEBUG
endif

sdssel: SDSSel.cpp 
	g++ -fopenmp -std=gnu++0x $(OPTFLAG) -Wall SDSSel.cpp -o sdssel
clean:
	rm -f sdssel
