OBJSD = sourceDew.o date_time.o date.o ctime.o 
OBJSP = sourcePreDew.o date_time.o date.o ctime.o
C++C=g++

dew: 	$(OBJSD)
	$(C++C) $(C++FLAGS) $(OBJSD) -o dew 

predew:	$(OBJSP)
	$(C++C) $(C++FLAGS) $(OBJSP) -o predew 

stationMask:	sourceStationMask.cpp 
		$(C++C) $(C++FLAGS) sourceStationMask.cpp -o stationMask

mask:	sourceMask.cpp 
		$(C++C) $(C++FLAGS) sourceMask.cpp -o mask

riverFlow:	sourceRiverFlow.cpp
		$(C++C) $(C++FLAGS) sourceRiverFlow.cpp -o riverFlow

clean:
	rm -f $(OBJSD) $(OBJSP) *.o core
