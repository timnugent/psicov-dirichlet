CFLAGS=-m64 -O3 -Iinclude
PSIFLAGS=-Wno-write-strings -mfpmath=sse -msse3 -funroll-loops
LIBS=-lz -lm -lgfortran

psicov:
	gfortran -c src/glasso_psicov.f90 -o src/glasso_psicov.o
	g++ $(CFLAGS) $(PSIFLAGS) -c src/psicov.cc -o src/psicov.o
	g++ $(CFLAGS) -c src/Alph.cc -o src/Alph.o
	g++ $(CFLAGS) -c src/Alphabet.cc -o src/Alphabet.o
	g++ $(CFLAGS) -c src/AlphabetTuple.cc -o src/AlphabetTuple.o
	g++ $(CFLAGS) -c src/DirichletReg.cc -o src/DirichletReg.o
	g++ $(CFLAGS) -c src/Filenames.cc -o src/Filenames.o
	g++ $(CFLAGS) -c src/gzstream.cc -o src/gzstream.o
	g++ $(CFLAGS) -c src/Input.cc -o src/Input.o
	g++ $(CFLAGS) -c src/LogGamma.cc -o src/LogGamma.o
	g++ $(CFLAGS) -c src/Multinomial.cc -o src/Multinomial.o
	g++ $(CFLAGS) -c src/NamedClass.cc -o src/NamedClass.o
	g++ $(CFLAGS) -c src/NamedObject.cc -o src/NamedObject.o
	g++ $(CFLAGS) -c src/NameToPtr.cc -o src/NameToPtr.o
	g++ $(CFLAGS) -c src/OrderSort.cc -o src/OrderSort.o
	g++ $(CFLAGS) -c src/Regularizer.cc -o src/Regularizer.o
	g++ $(CFLAGS) -c src/zfstream.cc -o src/zfstream.o
	g++ $(CFLAGS) -o psicov_d src/*.o $(LIBS)

clean:
	rm bin/psicov_d src/*.o