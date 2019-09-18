/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo and Clément Courbet
*
* This file is part of PPMC.
*
* PPMC is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* PPMC is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with PPMC.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/
#include <compression/ppmc.h>
#include <unistd.h>
using namespace CGAL;
using namespace std;

int main(int argc, char** argv) {

	if(argc<2){
		//here we only take the prefix of output2 as the parameter
		//for judging which dataset we are in
		std::cerr<<"usage: ppmc /path/to/second/dataset"<<std::endl;
		return 0;
	}

	/*firstly we need to get which data set are we processing*/
	char *prefix_2 = argv[1];
	int join_idx = -1; // index of the current file (0 or 1) matching to dataset 1 or 2
	char* mapper_id = getenv("mapreduce_task_id");
	char *stdin_file_name = getenv("mapreduce_map_input_file");// name of the input file

	/* This happens if program is not run in mapreduce
	 *  For testing locally, set/export the environment variable above */
	if (!stdin_file_name) {
		#ifdef DEBUG
		std::cerr << "Environment variable mapreduce_map_input_file is not set correctly." << std::endl;
		#endif
		return -1;
	}
	if(!mapper_id){
		#ifdef DEBUG
		std::cerr << "Environment variable mapreduce_task_id is not set correctly." << std::endl;
		#endif
		return -1;
	}
	#ifdef DEBUG
	std::cerr <<"stdin file name: "<< stdin_file_name<<"\nmapper_id:" <<mapper_id<< std::endl;
	#endif

	if (prefix_2!=NULL && strstr(stdin_file_name, prefix_2) != NULL) {
		join_idx = SID_2;
	} else {
		join_idx = SID_1;
	}

	if (join_idx < 0) {
		#ifdef DEBUG
		std::cerr << "Invalid join index" << std::endl;
		#endif
		return -1;
	}

	//process objects line by line and generate compressed data in
	//a temporary file and mbbs in hdfs
	std::stringstream output_path;
	output_path << "/tmp/compressionoutput" << mapper_id;
	#ifdef DEBUG
	std::cerr<<"compress data into "<<output_path.str()<<std::endl;
	#endif
	if (!compress_data(output_path.str(), mapper_id, join_idx)) {
		std::cerr << "Error reading input in" << std::endl;
		return -1;
	}
	#ifdef DEBUG
	Timer t;
	#endif

	return 0;
}


bool compress_data( std::string output_path, char* mapper_id, long join_id){

	MyMesh *currentMesh = NULL;
	long count_objects = -1;
	long obj_id = 0;
	std::string input_line; // Temporary line
	std::vector<std::string> fields; // Temporary fields
	std::stringstream ss;
	
	int i_mode = COMPRESSION_MODE_ID; // compression mode

	// Codec features status.
	bool b_useAdaptiveQuantization = false;
	bool b_useLiftingScheme = true;
	bool b_useCurvaturePrediction = true;
	bool b_useConnectivityPredictionFaces = true;
	bool b_useConnectivityPredictionEdges = true;
	bool b_allowConcaveFaces = true;
	bool b_useTriangleMeshConnectivityPredictionFaces = true;
	unsigned i_quantBit = 12;
	unsigned i_decompPercentage = 100;
	
	// MBB and space info
	double tmp_x, tmp_y, tmp_z;
	double low[3];
	double high[3];
	bool firstLineRead = false;
	/* Space info */
	double space_low[3];
	double space_high[3];


	// output binary file for the compressed data
	int offset = 0;
	//ofstream myFile (output_path, ios::out | ios::binary);
	int count = 0;
	char *buffer = new char[20000000];
	char *meshbuffer = new char[BUFFER_SIZE];
    for (size_t i = 0; i < BUFFER_SIZE; ++i) {
    	meshbuffer[i] = 0;
   	}
	char *currentPos = buffer;
	while (std::cin && getline(std::cin, input_line) && !std::cin.eof()) {
		// the input format is "ID	OFF|numbers|numbers|...."
		// return character in OFF file is replaced with bar "|"

		tokenize(input_line, fields, TAB, true);
		std::stringstream buffer(fields[0]);
		buffer >> obj_id;

		/* Parsing polyhedron input */
		try {
			// convert it back to a normal OFF format
			boost::replace_all(fields[1], BAR, "\n");
			// Init the random number generator.
			srand(4212);
			//read the mesh:
			//global variables
			currentMesh = new MyMesh(NULL, i_decompPercentage,
					     i_mode, i_quantBit, b_useAdaptiveQuantization,
					     b_useLiftingScheme, b_useCurvaturePrediction,
					     b_useConnectivityPredictionFaces, b_useConnectivityPredictionEdges,
					     b_allowConcaveFaces, b_useTriangleMeshConnectivityPredictionFaces,
					     fields[1], NULL, 0, meshbuffer);
			currentMesh->completeOperation();
			
			// output the mbb information
			// note that the first letter "0", same as the word "SPACE"
			// below, is the folder for this specific output
			std::cout << "0" << TAB <<  mapper_id << TAB << obj_id << TAB << join_id
					  << TAB << currentMesh->bbMin0.x() << TAB << currentMesh->bbMin0.y() << TAB << currentMesh->bbMin0.z()
					  << TAB << currentMesh->bbMax0.x() << TAB << currentMesh->bbMax0.y() << TAB << currentMesh->bbMax0.z()
					  << TAB << offset << TAB << currentMesh->dataOffset << std::endl; // offset is the beginning point of this object
			
			// write the compressed data into binary file
			memcpy(currentPos, currentMesh->p_data, currentMesh->dataOffset);
    		currentPos = currentPos + currentMesh->dataOffset;

			// compute offset for next object
			offset += currentMesh->dataOffset;

			// Get or update the space info
			low[0] = currentMesh->bbMin0.x();
			low[1] = currentMesh->bbMin0.y();
			low[2] = currentMesh->bbMin0.z();
			high[0] = currentMesh->bbMax0.x();
			high[1] = currentMesh->bbMax0.y();
			high[2] = currentMesh->bbMax0.z();
			if (!firstLineRead) {
				space_low[0] = low[0];
				space_low[1] = low[1];
				space_low[2] = low[2];
				space_high[0] = high[0];
				space_high[1] = high[1];
				space_high[2] = high[2];
				firstLineRead = true;
			} else {
				space_low[0] = low[0] < space_low[0] ? low[0] : space_low[0];
				space_low[1] = low[1] < space_low[1] ? low[1] : space_low[1];
				space_low[2] = low[2] < space_low[2] ? low[2] : space_low[2];
				space_high[0] = high[0] > space_high[0] ? high[0] : space_high[0];
				space_high[1] = high[1] > space_high[1] ? high[1] : space_high[1];
				space_high[2] = high[2] > space_high[2] ? high[2] : space_high[2];
			}
		} catch (...) {
			std::cerr << "******Geometry Parsing Error******" << std::endl;
			return -1;
		}
	
		fields.clear();
		delete currentMesh;
		count_objects++;
	} // end of while

	// output the overall space information
	// TODO move the space information into the binary file
	// and processed by the combine tool
	std::cout << "SPACE" << TAB << "T" << TAB
			  << space_low[0] << TAB << space_low[1] << TAB << space_low[2] << TAB
			  << space_high[0] << TAB << space_high[1] << TAB << space_high[2] << TAB
			  << count_objects << std::endl;
	
	// write the compressed data
	const char *cstr = output_path.c_str();
	std::ofstream myFile (cstr, std::ios::out | std::ios::binary);

	int amt = 67108864;
	long wridx = 0;
	for (wridx = 0; wridx + amt < offset; wridx = wridx + amt) {
		myFile.write((buffer + wridx), amt);
	}
	// writing last chunk
	if (offset > 0 && wridx < offset) {
		myFile.write( (buffer + wridx), offset % amt);
	}
	myFile.close();

	delete[] buffer;
    delete[] meshbuffer;

	#ifdef DEBUG
    std::cerr <<"processed "<<count_objects<<" objects"<<endl;
	std::cerr <<"total size of compressed data is "<<offset << std::endl; // the total size
	#endif

	return true;
}

