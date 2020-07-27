CFLAGS=`pkg-config --cflags opencv`
LDFLAGS=`pkg-config --libs-only-L opencv` -lopencv_core -lopencv_highgui -lopencv_imgcodecs -lopencv_imgproc
ssimulacra: ssimulacra.cpp
	g++ -std=c++11 -O2 -fstrict-aliasing -ffast-math $(CFLAGS) ssimulacra.cpp $(LDFLAGS) -o ssimulacra
clean:
	rm ssimulacra
install:
	cp ssimulacra /usr/local/bin/