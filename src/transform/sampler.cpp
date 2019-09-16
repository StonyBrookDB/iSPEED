#include <iostream>
#include <cstring>
#include <cmath>
#include <cstdlib> 
#include <limits>
#include <random>
#include <ctime>

using namespace std;

/*
 * This program performs random coin flipping with a fixed sampling rate
 * take input from standard in, and emit to standard out with a
 * sample rate
 * */

int main(int argc, char **argv) {

	double sample_rate = 1;
	if (argc >= 2) {
		sample_rate = atof(argv[1]);
	}

	if(sample_rate > 1 || sample_rate <= 0){
		std::cerr<<"sample rate should be a number between 0 and 1"<<std::endl;
	}

	int max = std::numeric_limits<int>::max();
	int threshold = static_cast<int>(ceil(max * sample_rate));

	srand(time(NULL));

	#ifdef DEBUG
	cerr << "sample rate: " << argv[1] << endl;
	#endif
	string input_line = "";
	int total = 0;
	int emitted = 0;
	while(cin && getline(cin, input_line) && !cin.eof()) {
		if (rand() % max <= threshold ) {
			cout << input_line << endl;
			emitted++;
		}
		total++;
	}
	cout.flush();
	#ifdef DEBUG
	cerr << emitted <<" out of "<<total << " is sampled"<< endl;
	#endif
}
