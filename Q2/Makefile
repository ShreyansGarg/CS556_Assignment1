all: q2 run-cpp run-python

q2: q2.cpp
	g++ q2.cpp -o q2

run-cpp: q2
	./q2

run-python:
	python3 -u plot-call-center-metrics.py

clean:
	rm -f q2
