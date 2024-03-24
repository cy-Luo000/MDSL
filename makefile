all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o KDtest main.cpp Graph.cpp -w
debug:
	g++ -g -O0 -march=native -o KDtest main.cpp Graph.cpp -w
sub:
	g++ -flto -Ofast -DNDEBUG -o MCC SubGraph.cpp
clean:
	rm -rf KPLEX
test:
	g++ -flto -Ofast -DNDEBUG -march=native -o KDtest main.cpp Graph.cpp -w
	./KDtest ../socfb/socfb-Berkeley13.bin 1
	./ZhouBD ../socfb/socfb-Berkeley13.bin 1
	./KDtest ../socfb/socfb-Texas84.bin 1
	./ZhouBD ../socfb/socfb-Texas84.bin 1
	./KDtest ./data/soc-youtube.bin 1
	./ZhouBD ./data/soc-youtube.bin 1