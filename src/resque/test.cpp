/*
 * test.cpp
 *
 *  Created on: Sep 25, 2019
 *      Author: teng
 */
#include "resque_3d.hpp"

#include <sys/time.h>
struct timeval get_cur_time(){
	struct timeval t1;
	gettimeofday(&t1, NULL);
	return t1;
}
double get_time_elapsed(struct timeval t1){
	struct timeval t2;
    double elapsedTime;
	gettimeofday(&t2, NULL);
	// compute and print the elapsed time in millisec
	elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
	elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
	return elapsedTime;
}

bool optimization = false;
using namespace std;
bool b_useAdaptiveQuantization = optimization;
bool b_useLiftingScheme = optimization;
bool b_useCurvaturePrediction = optimization;
bool b_useConnectivityPredictionFaces = optimization;
bool b_useConnectivityPredictionEdges = optimization;
bool b_allowConcaveFaces = true;
bool b_useTriangleMeshConnectivityPredictionFaces = true;
unsigned i_quantBit = 12;
unsigned i_decompPercentage = 100;

MyMesh *compress(string input){
	int i_mode = COMPRESSION_MODE_ID; // compression mode

	srand(PPMC_RANDOM_CONSTANT);
	MyMesh *currentMesh = new MyMesh(i_decompPercentage,
		 i_mode, i_quantBit, b_useAdaptiveQuantization,
		 b_useLiftingScheme, b_useCurvaturePrediction,
		 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
		 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
		 (char*)(input.c_str()), input.size());
	currentMesh->completeOperation();
	return currentMesh;
}

void print_mesh(MyMesh *mesh){
	std::stringstream os;
	os << *mesh;
	cout << os.str()<<endl;
}

void print_mesh_file(MyMesh *mesh, char *path){
	ofstream myfile;
	myfile.open(path);
	myfile << *mesh;
	myfile.close();
}

int main()
{
	string input_line;
	MyMesh *compressed = NULL;
	int index = 0;
	long total = 0;
	while(std::cin && getline(std::cin, input_line) && !std::cin.eof()){
		struct timeval t1 = get_cur_time();
		boost::replace_all(input_line, BAR, "\n");
		compressed = compress(input_line);
		total += compressed->dataOffset;
		index += 10;
		std::cerr<<"processed "<<index/10<<" "<<get_time_elapsed(t1)/1000<<endl;
		break;
	}
	char path[100];
	for(int i=0;i<=10;i++){
		srand(PPMC_RANDOM_CONSTANT);
		MyMesh *decompressed = new MyMesh(10*i,
				DECOMPRESSION_MODE_ID, i_quantBit, b_useAdaptiveQuantization,
						 b_useLiftingScheme, b_useCurvaturePrediction,
						 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
						 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
						 compressed->p_data, compressed->dataOffset);
		decompressed->completeOperation();
		sprintf(path,"lod%d.off",i*10);
		print_mesh_file(decompressed, path);

		cerr<<index<<" "<<decompressed->dataOffset<<endl;
		delete decompressed;

	}
//	int i_mode = DECOMPRESSION_MODE_ID; // compression mode
//	for(int i = 1;i<=10;i++){
//		int i_decompPercentage = 10*i;
//		srand(PPMC_RANDOM_CONSTANT);
//		MyMesh *decompressed = new MyMesh(i_decompPercentage,
//						 i_mode, i_quantBit, b_useAdaptiveQuantization,
//						 b_useLiftingScheme, b_useCurvaturePrediction,
//						 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
//						 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
//						 compressed->p_data, compressed->dataOffset);
//		decompressed->completeOperation();
//		cerr<<10*i<<" "<<decompressed->dataOffset<<endl;
//		delete decompressed;
//	}
//
//	for(int i=10;i<=100;i+=10){
//		struct timeval start = get_cur_time();
//		//repeat 10 times for profiling decompression performance
//		for(int r = 0;r<10;r++){
//			MyMesh *decompressed = new MyMesh(i,
//							 i_mode, i_quantBit, b_useAdaptiveQuantization,
//							 b_useLiftingScheme, b_useCurvaturePrediction,
//							 b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
//							 b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
//							 compressed->p_data, compressed->dataOffset);
//			decompressed->completeOperation();
//			delete decompressed;
//		}
//		cerr<<"decompression takes "<<get_time_elapsed(start)/10<<endl;
//	}
	delete compressed;

	return EXIT_SUCCESS;
}
