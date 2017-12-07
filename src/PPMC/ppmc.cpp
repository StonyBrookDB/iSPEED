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
#include <PPMC/ppmc.h>
#include <unistd.h>
using namespace CGAL;

//global variables
MyMesh *currentMesh = NULL;

// main method
// Initialize default values in spatial operator structure
void init(struct query_op &stop, struct query_temp &sttemp) {
	stop.extract_mbb = true;
	stop.collect_mbb_stat = false;
	stop.use_sampling = false;
	stop.sample_rate = 1.0;
	stop.offset = 0;

	stop.prefix_1 = NULL;
	stop.prefix_2 = NULL;
	stop.shape_idx_1 = 0;
	stop.shape_idx_2 = 0;
}

int main(int argc, char** argv) {

	struct query_op stop;
        struct query_temp sttemp;


        init(stop, sttemp);


	if (!extract_params(argc, argv, stop, sttemp)) {
		#ifdef DEBUG 
		std::cerr <<"ERROR: query parameter extraction error." << std::endl 
			<< "Please see documentations, or contact author." << std::endl;
		#endif
		usage();
		return -1;
	}

	char* stdin_file_name = NULL; // name of the input file
	int join_idx = -1; // index of the current file (0 or 1) matching to dataset 1 or 2
	int geom_idx = -1; // geometry field index

	//std::cout.precision(20);

	char* mapper_id = getenv("mapreduce_task_id");

	stdin_file_name = getenv("mapreduce_map_input_file");

	//std::cout << stdin_file_name << std::endl;

	if (!stdin_file_name) {
		/* This happens if program is not run in mapreduce
		 *  For testing locally, set/export the environment variable above */
		#ifdef DEBUG
		std::cerr << "Environment variable mapreduce_map_input_file is not set correctly." << std::endl;
		#endif
		return -1;
	}

	std::string input_line;
	/*while (cin && getline(cin, input_line) && !cin.eof()) {
    		std::cout << "MBB" << TAB << stdin_file_name << TAB << mapper_id << std::endl;
	}
	return 0;
	*/

	if (strstr(stdin_file_name, stop.prefix_1) != NULL) {
		join_idx = SID_1;
		geom_idx = stop.shape_idx_1;
	} else if (strstr(stdin_file_name, stop.prefix_2) != NULL) {
		join_idx = SID_2;
		geom_idx = stop.shape_idx_2;
	} else {
		std::cerr << "File name from environment variable \
			does not match any path prefix" << std::endl;
	}

	if (join_idx < 0) {
		#ifdef DEBUG
		std::cerr << "Invalid join index" << std::endl;
		#endif
		return -1;
	}

	//char* output_file_name = argv[1];
	//std::cout << output_file_name << std::endl;
	
	//std::string mapperidstr(mapper_id);
//	std::string output_path = "/tmp/compressionoutput" + mapperidstr; 
	std::stringstream output_path;
	output_path << "/tmp/compressionoutput" << mapper_id;
//	output_path << "/scratch/hadoopgis3d/mrbin/compressionoutput" << mapper_id;
			
	//std::string output_path = "/scratch/hadoopgis3d/mrbin/compressionoutput" + mapperidstr;  

	if (!compress_data(stdin_file_name, output_path.str(), mapper_id, join_idx)) {
		std::cerr << "Error reading input in" << std::endl;
		return -1;
	}
	#ifdef DEBUGTIME
		Timer t; 
	#endif

	return 0;
}


bool compress_data(char* stdin_file_name, std::string output_path, char* mapper_id, long join_id)
{
	long count_objects = -1;
	long obj_id = 0;
	std::string input_line; // Temporary line
	std::vector<std::string> fields; // Temporary fields
	std::stringstream ss;
	
	int i_mode = COMPRESSION_MODE_ID; // compression mode

	// Codec features status.
	bool b_useAdaptiveQuantization = false;
	//bool b_useAdaptiveQuantization = false;
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
		
		// Removal of \r symbol on Windows
		if (input_line.at(input_line.size() - 1) == '\r') {
			input_line = input_line.substr(0, input_line.size() - 1);
		}
		tokenize(input_line, fields, TAB, true);

		std::stringstream buffer(fields[0]);
		buffer >> obj_id;
		//obj_id = stod(fields[0]);

		count_objects++;
		//std::cout << count_objects << std::endl;
		/* Parsing polyhedron input */
		try {
			boost::replace_all(fields[1], BAR, "\n");
			//ss(fields[1]);	
			//std::cout << input_line << std::endl;
			//ss >> *poly;

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

			// Run the complete job and exit.
			//CGAL::Timer user_time;
			//user_time.start();  
			currentMesh->completeOperation();
			//std::cerr << "Execution       : " << user_time.time() << " seconds." << std::endl;
			
			
			//std::cout << "IDX"
			//	<< TAB 

			// for local test			
			std::cout << "0" << TAB <<  mapper_id << TAB << obj_id << TAB << join_id
				<< TAB << currentMesh->bbMin0.x() << TAB << currentMesh->bbMin0.y() << TAB << currentMesh->bbMin0.z() 
				<< TAB << currentMesh->bbMax0.x() << TAB << currentMesh->bbMax0.y() << TAB << currentMesh->bbMax0.z() 
				<< TAB << offset << TAB << currentMesh->dataOffset << std::endl; // offset is the beginning point of this object
			
			// write the comprssed data into binary file
			/*std::stringstream ss1;
			ss1 << "/scratch/hadoopgis3d/mrbin/compressionoutput" << count;
			count++;  
			ofstream myFile1 (ss1.str(), ios::out | ios::binary);*/
    			//myFile.write(currentMesh->p_data, currentMesh->dataOffset);
    			//
    			//
			//memcpy((buffer + offset), currentMesh->p_data, currentMesh->dataOffset);
			memcpy(currentPos, currentMesh->p_data, currentMesh->dataOffset);
    			currentPos = currentPos + currentMesh->dataOffset;
			//std::cerr << obj_id << TAB << ((const void *) buffer) << TAB << ((const void *) currentPos) << TAB << currentMesh->dataOffset << TAB << ((const void *) currentMesh->p_data) << std::endl;
			//sleep(3);
			//myFile1.close();
			//std::cout << "Content of first byte: " << (* ((int *)(currentMesh->p_data)) ) << TAB <<  (* ((int *)(currentMesh->p_data + sizeof(int))) ) << std::endl;

			// compute offset for next object
			offset += currentMesh->dataOffset;

			// Get the space info
			/* Collecting information about the space dimension */
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
		}
		catch (...) {
			std::cerr << "******Geometry Parsing Error******" << std::endl;
			return -1;
		}

		#ifdef DEBUG
		std::cerr << "Processing " << count_objects << std::endl;
		#endif
	
		fields.clear();
		delete currentMesh;
	} // end of while

	/* Output dimensions of space */
	std::cout << "SPACE" << TAB << "T" << TAB << space_low[0] << TAB << space_low[1] << TAB << space_low[2] 
			<< TAB << space_high[0]	<< TAB << space_high[1]	<< TAB << space_high[2] << TAB << count_objects << std::endl;
	
	const char *cstr = output_path.c_str();
	std::ofstream myFile (cstr, std::ios::out | std::ios::binary);
	//ofstream myFile (output_path, ios::out | ios::binary);

	
	// Writing in chunks
	int amt = 67108864;
	long wridx = 0;
	for (wridx = 0; wridx + amt < offset; wridx = wridx + amt) {
		myFile.write((buffer + wridx), amt);
	}
	// writing last chunk
	if (offset > 0 && wridx < offset) {
		myFile.write( (buffer + wridx), offset % amt);
	}
	
	// Old code: 
	// myFile.write(buffer, offset);
	myFile.close();

	delete[] buffer;
        delete[] meshbuffer;
	#ifdef DEBUG
	std::cerr << offset << std::endl; // the total size
	#endif
        std::cerr << "size of types: " << sizeof(uint2) << TAB << sizeof(uint4) << TAB << sizeof(uint) << std::endl;
	return true;

}

