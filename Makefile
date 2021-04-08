#############################################################################################################
#
# MAKEFILE for Approximated-FQN 
#
# @brief Online Qn estimator using a uniformly collapsed sketch
#
# MODE
#
# -DTEST used to perform only processing
# -DCHECK used to log exact and estimated quantile and median, along with outliers and inliers
#
# @note 
# -DCMP to log all metrics into a single file for comparison against FQN
#############################################################################################################


OS=$(shell uname -s)
ifeq ($(OS),Linux)
	CC=icc
	CFLAGS=-std=c++14 -O3 -DCMP
else
	CC=clang++
	CFLAGS=-std=c++14 -Os -DCMP
endif


TARGET=AFQN7
DEPS=src/IIS.cc src/QuickSelect.cc src/Utility.cc src/DDSketch.cc src/Approx-FQN-Test.cc

LDFLAGS=


MODE=-DTEST#-DCHECK #


all:$(TARGET)

$(TARGET):
	@echo "Compiling for " $(OS)
	$(CC) $(CFLAGS) -o $(TARGET) $(DEPS) $(MODE) $(DIFFS) $(SAMPLE) $(LDFLAGS)



clean:
	rm -f *~ $(TARGET) log.txt err.txt *.csv
	rm -rf $(TARGET).dSYM
	