target   := ./rec
obj_dir  := ./objects/
sources  := $(shell find ./sources -name \*.cpp)
objects  := $(addprefix $(obj_dir), $(notdir $(sources:.cpp=.o)))
dirs     := $(dir $(sources))
VPATH    := $(dirs)

cxx = g++
cxxflags = -g -O2 -std=c++11 -I/usr/local/include/eigen3
lflags = -std=c++17

$(target): $(objects)
	$(cxx) -o $@ $^ $(lflags)

$(objects): $(obj_dir)%.o: %.cpp
	$(cxx) -o $@ -c $< $(cxxflags)

clean:
	rm -f $(target)
	rm -f $(obj_dir)*.o
