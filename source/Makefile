TARGETS=support
LIBS=readobj.cpp overhang.cpp

dummy: $(TARGETS)
	echo end

% : %.cpp
	g++ $< -o $@ $(LIBS)