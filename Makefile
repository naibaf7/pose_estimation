CXX = g++
CXXFLAGS = -Wall -std=c++11

UNAME := $(shell uname -s)
ifeq ($(UNAME), Linux)
	LINUX := 1
else ifeq ($(UNAME), Darwin)
	OSX := 1
endif

ifneq ($(OSX), 1)
	CXXFLAGS += -fopenmp
endif

CXXRUN = -O3
CXXDBG = -O0 -g -DDBG

FILES = P3p.cpp P4pf.cpp parse_bundler.cpp bundler_camera.cpp SIFT_loader.cpp pose_utils.cpp query_loader.cpp query_processor_basic.cpp query_processor_advanced.cpp query_processor.cpp import_export.cpp benchmark.cpp

INCLUDE = 
LIBRARY = -lopencv_core -lopencv_highgui -lopencv_imgproc

all : pose_estimation pose_estimation_dbg p3p_testing p4p_testing

pose_estimation : pose_estimation.cpp $(FILES)
	$(CXX) $(CXXFLAGS) $(CXXRUN) $(INCLUDE) -o pose_estimation pose_estimation.cpp $(FILES) $(LIBRARY) 	

pose_estimation_dbg : pose_estimation.cpp $(FILES)
	$(CXX) $(CXXFLAGS) $(CXXDBG) $(INCLUDE) -o pose_estimation_dbg pose_estimation.cpp $(FILES) $(LIBRARY)

p3p_testing: p3p_testing.cpp P3p.cpp P4pf.cpp
	$(CXX) $(CXXFLAGS) $(CXXDBG) $(INCLUDE) -o p3p_testing p3p_testing.cpp P3p.cpp P4pf.cpp $(LIBRARY)

p4p_testing: p4p_testing.cpp P4pf.cpp
	$(CXX) $(CXXFLAGS) $(CXXDBG) $(INCLUDE) -o p4p_testing p4p_testing.cpp P4pf.cpp $(LIBRARY)

clean:
	rm -r pose_estimation pose_estimation_dbg p3p_testing