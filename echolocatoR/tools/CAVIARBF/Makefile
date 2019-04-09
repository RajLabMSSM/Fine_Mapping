all: model_search caviarbf
model_search: src/model_search.cpp src/bf_io.cpp src/InputArgument.cpp src/CombinationTable.cpp src/model_selection.cpp
	cd ./src; g++ -Wall -Wno-comment -O3 -I../tclap-1.2.1/include model_search.cpp bf_io.cpp InputArgument.cpp CombinationTable.cpp model_selection.cpp -o ../model_search

caviarbf: src/InputArgument.cpp src/calculateBayesFactors.cpp src/caviarbf.cpp 
	cd ./src; g++ -Wall -Wno-comment -O3 -I../eigen-eigen-1306d75b4a21/ -I../tclap-1.2.1/include/ InputArgument.cpp calculateBayesFactors.cpp caviarbf.cpp -o ../caviarbf
