all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o QDC-init main.cpp Graph.cpp -w
debug:
	g++ -g -O0 -march=native -o QDC-init main.cpp Graph.cpp -w
sub:
	g++ -flto -Ofast -DNDEBUG -o MCC SubGraph.cppQ
clean:
	rm -rf QDC-init
test:
	g++ -flto -Ofast -DNDEBUG -march=native -o QDC-init main.cpp Graph.cpp -w
	./KDC-init ../socfb/socfb-Berkeley13.bin 1
	./ZhouBD ../socfb/socfb-Berkeley13.bin 1
	./KDC-init ../socfb/socfb-Texas84.bin 1
	./ZhouBD ../socfb/socfb-Texas84.bin 1
	./KDC-init ./data/soc-youtube.bin 1
	./ZhouBD ./data/soc-youtube.bin 1