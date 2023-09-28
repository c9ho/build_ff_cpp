CXX = g++

CXXFLAGS = -std=c++2a

SRCS = build_ff_main.cpp ForceField.cpp ForceFieldAngle.cpp ForceFieldAtom.cpp ForceFieldBond.cpp ForceFieldDihedral.cpp ForceFieldImproper.cpp ForceFieldVDW.cpp functions.cpp Guest.cpp MetalPairsReader.cpp Processor.cpp ReferenceConfig.cpp tool.cpp

OBJS = $(SRCS:.cpp=.o)

TARGET = build_force_field

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)