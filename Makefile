CC = clang++
CFLAGS = -g -Wall -Ofast
LIBS = -DEIGEN_DONT_PARALLELIZE -fopenmp=libomp
OPTS = -fsave-optimization-record -foptimization-record-file=./opt-viewer/opts.yaml

TARGET = local_search

.PHONY: $(TARGET)
all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(LIBS) $(TARGET).cpp -o $(TARGET)

clean:
	$(RM) $(TARGET)

.PHONY: opt
opt:
	cp *.h ./opt-viewer
	cp *.cpp ./opt-viewer
	python3 ./opt-viewer/opt-viewer.py ./opt-viewer/opts.yaml -o ./opt-viewer
	
