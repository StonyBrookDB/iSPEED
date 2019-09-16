/* 
 * This is the framework controller program. It manages calls to MapReduce and holds results, parameters.
 * This serves as the basic front-end interface
 * */

#include <framework/queryprocessor_3d.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	cout.precision(20);
	struct framework_vars fr_vars;

	
	/* Initialize default values */
	fr_vars.numreducers = 1;

	init_params(fr_vars);

	if (!extract_params(argc, argv, fr_vars)) {
		#ifdef DEBUG
		cerr << "Not valid arguments" << endl;
		#endif
		return 1;
	}
	
	#ifdef DEBUGTIME
	time_t start_exec_time, end_exec_time;
	double total_exec_time;
	time(&start_exec_time);
	#endif


	/////////////////////////////////////////////
	/* Process query */
	/////////////////////////////////////////////

	/////////////////////
	//  Containment ////
	///////////////////
	//////////////////////////////////////////
	//  Spatial join and nearest neighbor ////
	/////////////////////////////////////////
	if (fr_vars.query_type.compare(QUERYPROC_JOIN) == 0) {
		execute_spjoin(fr_vars);		
	}
	#ifdef DEBUGTIME
		time(&end_exec_time);
		total_exec_time = difftime(end_exec_time,start_exec_time);
		cerr << "********************************************" << endl;
		cerr << "Total execution time: " 
			<< total_exec_time
			<< " seconds." << endl;
		cerr << "********************************************" << endl;
	#endif

}

#ifdef DEBUGSTAT
// Handle tile statistics
void post_process_stat(char *tmpFile, stringstream &output) {
	string input_line;
	vector<string> fields;
	ifstream inputfile(tmpFile);
	stringstream ss;

	double low[2];
	double high[2];

	int tile_count = 0;
	long obj_count = 0;
	double margins = 0;

	long min_count = LONG_MAX;
	long max_count = 0;

	accumulator_set<double, stats<tag::variance> > acc;
	while(getline(inputfile, input_line)) {
		tokenize(input_line, fields, TAB, true);
		low[0] = stod(fields[1]);
		low[1] = stod(fields[2]);
		high[0] = stod(fields[3]);
		high[1] = stod(fields[4]);
		obj_count = strtol(fields[5].c_str(), NULL, 10);
		margins += (high[0] - low[0]) * (high[1] - low[1]);

		acc(static_cast<int>(obj_count));
		if (max_count < obj_count) {
			max_count = obj_count;
		}
		if (min_count > obj_count) {
			min_count = obj_count;
		}

		fields.clear();
		tile_count++;
	}
	output << tile_count << TAB
	//	<< spinfo.num_objects << TAB
		<< margins << TAB
		<< mean(acc) << TAB
		<< min_count << TAB
		<< max_count << TAB
		<< sqrt(variance(acc));
}
#endif

bool compress_data(string programpath, vector<string> &input_paths,
	string output_path, struct framework_vars &fr_vars) {
	hdfs_delete(programpath, output_path);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.task.timeout=36000000");

	arr_args.push_back("-libjars");
	arr_args.push_back(fr_vars.hadoopgis_prefix + CUSTOM_JAR_REL_PATH);
	arr_args.push_back("-outputformat");
	arr_args.push_back("com.custom.CustomMultiOutputFormat");
	
	/*
	arr_args.push_back("-input");	
	arr_args.push_back(input_path);
	*/
	for(vector<string>::iterator it = input_paths.begin() ; it != input_paths.end(); ++it) {
		arr_args.push_back("-input");
		arr_args.push_back(*it);
	}

	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + COMPRESSION);


	// should also work: 
	arr_args.push_back("-mapper");
	stringstream ss;
	ss << COMPRESSION << " " << fr_vars.input_path_2;
	arr_args.push_back(ss.str());
	ss.str("");
        //arr_args.push_back("cat");
	
	/*arr_args.push_back("-reducer");
	arr_args.push_back(SPACE_EXTRACTOR);

	// 2 reducers = 1 for outputting mbb and 1 for outputting spacial_dimension
	arr_args.push_back("-numReduceTasks");
	//arr_args.push_back("2");
	arr_args.push_back(fr_vars.numreducers_str);*/

	//arr_args.push_back("-numMapTasks");
	//
	/*
	int a = max(1, static_cast<int>(ceil(fr_vars.obtain_size_2/(5*1024*1024.0))));
	int b = fr_vars.numreducers;
	int c = max(1, static_cast<int>(ceil(fr_vars.obtain_size_2/(64*1024*1024.0))));
	int optimal = 1;
	if (a >= b && b >= c) {
		// everything else
		optimal = b;
	} else if (b >= a) {
		// extremely small data set
		optimal = a;
	} else {
		// large to very large data set
		optimal = c;
	}
	optimal = a; // 5M for the skeleton program 
	//arr_args.push_back(ss.str());
	*/
	arr_args.push_back("-numReduceTasks");
	arr_args.push_back("0");
	/*
	//arr_args.push_back("-jobconf");
	ss << "mapreduce.job.maps=" << optimal;
	arr_args.push_back(ss.str());
	*/
	//arr_args.push_back("-jobconf");

#ifdef DEBUG
	cerr << "Compress data program params: " << endl;
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(programpath, arr_args))) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in compressing data: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in compressing data: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}


bool partition_data(string programpath, string input_path,
	string output_path, string partitionmethod, int bucket_size,
	int step, double samplerate, struct framework_vars &fr_vars,
	char *cachefilefullpath) {
	char *cachefilename = strrchr(cachefilefullpath, '/'); // pointing to just the name of the cache file
	cachefilename++;
	hdfs_delete(programpath, output_path);
	stringstream ss;
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};

	arr_args.push_back("-input");
	arr_args.push_back(input_path);
	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	if (step == 1) {
		// For first step in partitioning
		arr_args.push_back(fr_vars.hadoopgis_prefix + MBB_SAMPLER);
	} else {
		// For 2nd step in partitioning
		arr_args.push_back(fr_vars.hadoopgis_prefix + MANIPULATE);
	}
	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + partitionmethod);

	if (cachefilename != NULL) {
		string strtmp(cachefilefullpath);
		arr_args.push_back("-file");
		arr_args.push_back(strtmp);
	}

	arr_args.push_back("-mapper");
	ss.str("");
	if (step == 1) {
		ss << MBB_SAMPLER << " " << samplerate;
	} else {
		ss << MANIPULATE << " " << cachefilename ;
		// no more mbb only mode	<< " -c " << cachefilename << " -m";
	}
	arr_args.push_back(ss.str());

	arr_args.push_back("-reducer");
	ss.str("");
	bucket_size = max(static_cast<int>(floor(bucket_size * samplerate)), 1);
	ss << partitionmethod << " -b " << bucket_size;
	if (cachefilename != NULL) {
		ss << " -c " << cachefilename;
	}
	arr_args.push_back(ss.str());

	arr_args.push_back("-numReduceTasks");
	if (step == 1) {
		arr_args.push_back("1");
	} else {
		//arr_args.push_back("0");
		arr_args.push_back(fr_vars.numreducers_str);
	}
	//arr_args.push_back(numreducers_str);

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");

#ifdef DEBUG
	cerr << "Partitioning params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(programpath, arr_args))) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in partitioning " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in partitioning: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}

// Read and combine all space info files
void read_space(char *filename, struct framework_vars &fr_vars) {
	string input_line;
	vector<string> fields;
	int pos = 0;
	bool firstLine = true;
	ifstream file(filename);
	int offset = 1;
	fr_vars.spinfo.num_objects = 0;
	while(getline(file, input_line)) {
		tokenize(input_line, fields, TAB, true);
		if (!firstLine) {
			fr_vars.spinfo.space_low[0] = min(fr_vars.spinfo.space_low[0], atof(fields[offset].c_str()));
			fr_vars.spinfo.space_low[1] = min(fr_vars.spinfo.space_low[1], atof(fields[offset + 1].c_str()));
			fr_vars.spinfo.space_low[2] = min(fr_vars.spinfo.space_low[2], atof(fields[offset + 2].c_str()));
			fr_vars.spinfo.space_high[0] = max(fr_vars.spinfo.space_high[0], atof(fields[offset + 3].c_str()));
			fr_vars.spinfo.space_high[1] = max(fr_vars.spinfo.space_high[1], atof(fields[offset + 4].c_str()));
			fr_vars.spinfo.space_high[2] = max(fr_vars.spinfo.space_high[2], atof(fields[offset + 5].c_str()));
			fr_vars.spinfo.num_objects += atol(fields[offset + 6].c_str());
		} else {
			// First line
			firstLine = false;
			fr_vars.spinfo.space_low[0] = atof(fields[offset].c_str());
			fr_vars.spinfo.space_low[1] = atof(fields[offset + 1].c_str());
			fr_vars.spinfo.space_low[2] = atof(fields[offset + 2].c_str());
			fr_vars.spinfo.space_high[0] = atof(fields[offset + 3].c_str());
			fr_vars.spinfo.space_high[1] = atof(fields[offset + 4].c_str());
			fr_vars.spinfo.space_high[2] = atof(fields[offset + 5].c_str());
			fr_vars.spinfo.num_objects += atol(fields[offset + 6].c_str());
		}	
		fields.clear();
	}
	file.close();
}



bool sp_join(string programpath, vector<string> &input_paths,
		string output_path, struct framework_vars &fr_vars,
		char *cachefilefullpath) {
	char *cachefilename = strrchr(cachefilefullpath, '/'); // pointing to just the name of the cache file
	cachefilename++; // Advance the pointer past the slash delimiter character

	hdfs_delete(programpath, output_path);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	//arr_args.push_back("-conf");
	//arr_args.push_back("/home/vhoang/mapred-site.xml");
	for(vector<string>::iterator it = input_paths.begin(); it != input_paths.end(); ++it) {
		arr_args.push_back("-input");
		arr_args.push_back(*it);
	}

	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + MANIPULATE);
	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + RESQUE);
	arr_args.push_back("-file");
	string strtmp(cachefilefullpath);
	arr_args.push_back(strtmp);

	// the mapper phase assign each object to
	// different tiles with its mbb
	arr_args.push_back("-mapper");
	stringstream ss;
	ss << MANIPULATE << " " << cachefilename;
	arr_args.push_back(ss.str());

	arr_args.push_back("-reducer");
	ss.str("");
	ss << RESQUE << " --lod " <<  fr_vars.decomp_lod<<" -j "<<fr_vars.join_cardinality;
	arr_args.push_back(ss.str()); // Offset to account for tile id and join index
	//arr_args.push_back("cat");

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back(fr_vars.numreducers_str);

	// mapper-only job
	//arr_args.push_back("0");

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");

#ifdef DEBUG
	cerr << "Executing spjoin program params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(programpath, arr_args))) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in sp join: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in sp join: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}


bool execute_spjoin(struct framework_vars &fr_vars) {
	// First check if dataset 1 has been loaded or not
	// Find the size of datasets 1 and 2
	//   in order to compute partitioning parameter (bucket size)
	// loaded_1 = hdfs_check_data(hadoopcmdpath, input_path_1 + "/" + PARTITION_FILE_NAME);
	fr_vars.obtain_size_1 = hdfs_get_size(fr_vars.hadoopcmdpath, fr_vars.input_path_1);
	if (fr_vars.obtain_size_1 >= 0) {
		fr_vars.spinfo.total_size = fr_vars.obtain_size_1;
	}

	// If there is a 2nd dataset, also check if it has been loaded and obtain its size
	if (fr_vars.join_cardinality > 1) {
		// loaded_2 = hdfs_check_data(hadoopcmdpath, input_path_2 + "/" + PARTITION_FILE_NAME);
		fr_vars.obtain_size_2 = hdfs_get_size(fr_vars.hadoopcmdpath, fr_vars.input_path_2);
		if (fr_vars.obtain_size_2 >= 0) {
			fr_vars.spinfo.total_size += fr_vars.obtain_size_2;
		}
	}
#ifdef DEBUG
	cerr << "Total size of 1 " << fr_vars.obtain_size_1 << endl;
	if (fr_vars.join_cardinality > 1) {
		cerr << "Total size of 2 " << fr_vars.obtain_size_2 << endl;
	}
	cerr << "Total object size: " << fr_vars.spinfo.total_size << endl;
	cerr << "Data 1 is loaded " << (fr_vars.loaded_1 ? "true" : "false" ) << endl;
	cerr << "Data 2 is loaded " << (fr_vars.loaded_2 ? "true" : "false" ) << endl;
#endif

	vector<string> inputpaths;

	/*// Neither data has been indexed */
	string mbb_output = fr_vars.output_path + "_mbb";
	inputpaths.push_back(fr_vars.input_path_1);
	if (fr_vars.join_cardinality > 1) {
		inputpaths.push_back(fr_vars.input_path_2);
	}
	// compress the input data
	// -- Extract object MBB and grab space dimension
	//  -- and compress data at the same time ([Hoang] thinks)
	cerr << "# of input paths: " << inputpaths.size() << TAB << fr_vars.join_cardinality << endl;
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.mbb_output)) {
		#ifdef DEBUG
		cerr << "\n\nCompressing & Extracting MBBs\n" << endl;
		#endif
		if (!compress_data(fr_vars.hadoopcmdpath, inputpaths,
				fr_vars.mbb_output, fr_vars)) {
			cerr << "Failed extracting MBB"  << endl;
			exit(1);
		}
	}

	//// Retrieve the total space dimension
	int tmpfd = mkstemp(nametemplate);
	char *tmpFile = nametemplate;
	close(tmpfd);
	tmpfd = open(tmpFile, O_RDWR | O_CREAT | O_TRUNC , 0777);
	#ifdef DEBUG
	cerr << "Temp file: " << tmpFile << endl;
	#endif

	// Create another file where you will write the content of all SPACE info paths into
	// Getting min_x, min_y, min_z, max_z, max_y, max_z
	// Obtain the cache file from hdfs
	bool res_cat = hdfs_cat(fr_vars.hadoopcmdpath, fr_vars.space_path, tmpfd);
	close(tmpfd);

	read_space(tmpFile, fr_vars);
	#ifdef DEBUG
	cerr << "Space dimensions: " << fr_vars.spinfo.space_low[0] << TAB
	<< fr_vars.spinfo.space_low[1] << TAB << fr_vars.spinfo.space_low[2] << TAB
	<< fr_vars.spinfo.space_high[0] << TAB << fr_vars.spinfo.space_high[1]
	<< TAB << fr_vars.spinfo.space_high[2] << endl;
	cerr << "Number objects: " << fr_vars.spinfo.num_objects << endl;
	#endif
	std::ofstream ofs;
	ofs.open (tmpFile, std::ofstream::out | std::ofstream::trunc);
	ofs << "T" << TAB << fr_vars.spinfo.space_low[0] << TAB << fr_vars.spinfo.space_low[1] << TAB << fr_vars.spinfo.space_low[2]
	 << TAB	<< fr_vars.spinfo.space_high[0] << TAB << fr_vars.spinfo.space_high[1] << TAB << fr_vars.spinfo.space_high[2]
		<< TAB << fr_vars.spinfo.num_objects << endl;
	ofs.close();

	// Saving those info into a struct that gets passed along query parameters into sequential steps


	if (fr_vars.comp_mode) {
		// stop the program
		cerr << "Done with compression" << endl;
		return 0;
	}

	// Run combiner on primary node and then loaders on all nodes
	//
	// runcombiner.sh kind of belongs here
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	//
	// Dejun:
	string allmbbspath = fr_vars.output_path + "_inputresque";
	vector<string> inputresque;
	inputresque.push_back(allmbbspath);

	cerr << "Start partitioning" << endl;

	// Setting default block size if it has not been set
	if (fr_vars.bucket_size < 0) {
		// Bucket size was not set
		double blockSize = 16000000; // approximately 16MB
		fr_vars.bucket_size = max(static_cast<int>(floor(blockSize
						/ fr_vars.spinfo.total_size * fr_vars.spinfo.num_objects)), 1);
	}

	// 1st step of partition the data to generate tile boundaries
	#ifdef DEBUG
	cerr << "Bucket size: " << fr_vars.bucket_size << endl;
	cerr << "Sampling rate: " << fr_vars.sampling_rate << endl;
	#endif

	if (!fr_vars.para_partition) {
		fr_vars.rough_bucket_size = fr_vars.bucket_size;
		if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.partitionpath)) {
			#ifdef DEBUG
			cerr << "\nPartitioning 1st steps\n" << endl;
			#endif
			// MapReduce partitioning happens here
			if (!partition_data(fr_vars.hadoopcmdpath, allmbbspath,
						fr_vars.partitionpath, fr_vars.partition_method, fr_vars.bucket_size,
						1, fr_vars.sampling_rate, fr_vars, tmpFile)) {
				cerr << "Failed partitioning 1st step" << endl;
				// Remove cache file
				remove(tmpFile);
				exit(1);
			}
		}
	} else {
		if (fr_vars.rough_bucket_size == -1) {
			// Rough bucket size has not been set
			if (fr_vars.numreducers > 0) {
				fr_vars.rough_bucket_size = max(static_cast<int>(floor(fr_vars.spinfo.num_objects
								/ fr_vars.numreducers / 2)), 1);
			} else {
				fr_vars.rough_bucket_size = max(static_cast<int>(floor(fr_vars.spinfo.num_objects
								/ 2)), 1);
			}
		}
		cerr << "Rough bucket size: " << fr_vars.rough_bucket_size << endl;
		if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.partitionpath)) {
			#ifdef DEBUG
			cerr << "\n\nPartitioning 1st step\n" << endl;
			#endif
			if (!partition_data(fr_vars.hadoopcmdpath, fr_vars.mbb_path,
						fr_vars.partitionpath, fr_vars.partition_method, fr_vars.rough_bucket_size,
						1, fr_vars.sampling_rate, fr_vars, tmpFile)) {
				cerr << "Failed partitioning 1st step" << endl;
				// Remove cache file
				remove(tmpFile);
				exit(1);
			}
		}
	}
	// Obtain the cache file from hdfs
	tmpfd = open(tmpFile, O_RDWR | O_CREAT | O_TRUNC , 0777);
	bool res_partition = hdfs_cat(fr_vars.hadoopcmdpath, fr_vars.partitionpathout, tmpfd);
	close(tmpfd);

	#ifdef DEBUG
	cerr << "Temp file name to hold partition boundary: " << tmpFile << endl;
	#endif

	if (fr_vars.para_partition) {
		// Second round of partitioning
		#ifdef DEBUG
		cerr << "\n\nPartitioning 2nd steps\n" << endl;
		#endif
		if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.partitionpath2)) {
			if (!partition_data(fr_vars.hadoopcmdpath, fr_vars.mbb_path,
						fr_vars.partitionpath2, fr_vars.partition_method_2, fr_vars.bucket_size,
						2, 1, fr_vars, tmpFile)) {
				cerr << "Failed partitioning 2nd step" << endl;
				// Remove cache file
				remove(tmpFile);
				exit(1);
			}
		}
		// Update/Overwrite the tmp file
		tmpfd = open(tmpFile, O_RDWR | O_CREAT | O_TRUNC , 0777);
		bool res_partition_2 = hdfs_cat(fr_vars.hadoopcmdpath, fr_vars.partitionpathout2, tmpfd);
		close(tmpfd);
	}

	// tmpFile contains partition index
#ifdef DEBUGSTAT
	//string stat_path = output_path + "_stat";
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, stat_path)) {
		if (!collect_stat(fr_vars.hadoopcmdpath, fr_vars.mbb_path, fr_vars.stat_path,
					fr_vars.sharedparams, tmpnameonly, tmpFile)) {
			cerr << "Failed obtaining stats" << endl;
			// Remove cache file
			remove(tmpFile);
			exit(1);
		}
	}
	cerr << "Done collecting tile counts" << endl;

	// Overwrite the cache file
	tmpfd = open(tmpFile, O_RDWR | O_CREAT | O_TRUNC , 0777);
	hdfs_cat(fr_vars.hadoopcmdpath, fr_vars.statpathout, tmpfd);
	close(tmpfd);

	stringstream outputss;
	outputss << (fr_vars.para_partition ? "true" : "false") << TAB
		<< fr_vars.partition_method << TAB
		<< fr_vars.rough_bucket_size << TAB
		<< fr_vars.sampling_rate << TAB
		<< fr_vars.partition_method_2 << TAB
		<< fr_vars.bucket_size << TAB
		<< fr_vars.spinfo.num_objects << TAB;
	post_process_stat(tmpFile, outputss);
	cout << outputss.str() << endl;

#else



	// Most important-heavy-lifting work here
	// inputresque: is the list of all objects MBBs with offset and length of the level 0 compression
	//  e.g. 1203 231 242341 12424 12441 42414 0 200 (starts at 0 byte and 200 bytes)
	// Spatial join step
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.joinoutputpath)) {
		#ifdef DEBUG
		cerr << "\n\nExecuting spatial joins\n" << endl;
		#endif
		if (!sp_join(fr_vars.hadoopcmdpath, inputresque, fr_vars.joinoutputpath, fr_vars,
					tmpFile)) {
			cerr << "Failed spatial join" << endl;
			// Remove cache file
			//
			//  Keep it for now
			// remove(tmpFile);
			exit(1);
		}
	}
	cerr << "Done with spatial join." << endl;

	// Perform duplicate removal

	/*	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.output_path)) {
#ifdef DEBUG
cerr << "\n\nBoundary object handling\n" << endl;
#endif
if (!duplicate_removal(fr_vars.hadoopcmdpath,
fr_vars.joinoutputpath, fr_vars.output_path, fr_vars)) {
cerr << "Failed boundary handling" << endl;
	// Remove cache file
	remove(tmpFile);
	exit(1);
	}
	}
	cerr << "Done with boundary handling. Results are stored at "
	<< fr_vars.output_path << endl;
	*/
#endif
	/*		cerr << "Cleaning up/Removing temporary directories" << endl;
			remove(tmpFile);
			if (fr_vars.remove_tmp_dirs) {
			hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.partitionpath);
			if (fr_vars.para_partition) {
			hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.partitionpath2);
			}
			hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.joinoutputpath);
			}
			if (fr_vars.remove_tmp_mbb) {
			hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.mbb_output);

			}
			*/

	cout.flush();
	cerr.flush();

	return true;
}




