#include "init.h"
#include "root.h"
#include "utils.h"

using namespace std;

#define TIMING time(&cur_time); cout << "TIME from the program start: " << cur_time - start_time << " sec\n";

//int main() {
int main(int argc, const char* argv[]) {
	
	time_t start_time, cur_time;
	time(&start_time);
	TIMING;
	
	if (doMainAnalysis) {
        analysis a;
        a.Loop();
	}
    if (0) {
        analysis a;
        a.LoopEle();
	}
    if (1 && doMainAnalysisPhoton) {
        analysis a;
        a.LoopPho();
	}
	if (doPostInitialResultsGraphs) {
		postInitialResultsAnalysis();
	}
	
	
	
} // END: main(int argc, const char* argv[]) {


