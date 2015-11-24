CC=gcc
CFLAGS=-c -Wall -I/usr/include/flycapture -I/usr/include/ -I/usr/local/include/opencv2 -I/usr/include/glibmm-2.4 -I/usr/lib64/glibmm-2.4/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include/ -I/usr/include/sigc++-2.0 -I/usr/lib64/sigc++-2.0/include/ -I/usr/include/boost 
LDFLAGS=   -L/usr/local/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_videoio -lsigc-2.0 -lglibmm-2.4 -lglib-2.0  -lstdc++ -lncurses -lflycapture
SOURCES=FrameRateCounter.cpp SaveVideoClass.cpp  SaveVideo.cpp 
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=savevideo


all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(DEBUG) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm FrameRateCounter.o SaveVideoClass.o SaveVideo.o savevideo
