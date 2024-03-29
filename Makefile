CXX       = g++
CXXFLAGS  = -g -O2 -Wall -c
LDFLAGS   = -g -O2 -Wall -L/usr/X11R6/lib -lGL -lglut
RM        = rm -f
EXEC      = minigl

all: $(EXEC)

$(EXEC): minigl.o main.o
	$(CXX) $(LDFLAGS) $^ -o $@

minigl.o: minigl.cpp minigl.h
	$(CXX) $(CXXFLAGS) $< -o $@

main.o: main.cpp minigl.h
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) *~ *.bak *.o $(EXEC)

onlocal: minigl.o
	$(CXX) -g -O2 -Wall main.cpp -L/usr/X11R6/lib -lGL -lglut minigl.o -o minigl
