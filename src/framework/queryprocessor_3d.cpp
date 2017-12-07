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
	fr_vars.numreducers = 10;

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

bool extract_mbb(string programpath, vector<string> &input_paths,
	string output_path, string original_params, struct framework_vars &fr_vars) {
	hdfs_delete(programpath, output_path);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};

	arr_args.push_back("-libjars");
	arr_args.push_back(fr_vars.hadoopgis_prefix + CUSTOM_JAR_REL_PATH);
	arr_args.push_back("-outputformat");
	arr_args.push_back("com.custom.CustomMultiOutputFormat");

	for(vector<string>::iterator it = input_paths.begin(); it != input_paths.end(); ++it) {
		arr_args.push_back("-input");
		arr_args.push_back(*it);
	}

	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + MANIPULATE);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + SPACE_EXTRACTOR);

	arr_args.push_back("-mapper");
	arr_args.push_back(MANIPULATE + " --offset 0" + original_params + " --extract ");

	arr_args.push_back("-reducer");
	arr_args.push_back(SPACE_EXTRACTOR);

	// 2 reducers = 1 for outputting mbb and 1 for outputting spacial_dimension
	arr_args.push_back("-numReduceTasks");
	//arr_args.push_back("2");
	arr_args.push_back(fr_vars.numreducers_str);

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");

	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

#ifdef DEBUG
	cerr << "Extract MBB program params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in extracting MBBs: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in extracting MBBs: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}

bool extract_skeleton(string programpath, string input_path,
	string output_path, struct framework_vars &fr_vars) {
	hdfs_delete(programpath, output_path);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.task.timeout=36000000");

	/*	
	arr_args.push_back("-Dmapreduce.max.split.size=10485760");
	arr_args.push_back("-Dmapred.max.split.size=10485760");
	arr_args.push_back("-Dmapreduce.min.split.size=10485760");
	arr_args.push_back("-Dmapred.min.split.size=10485760");
	arr_args.push_back("-Dmapreduce.input.fileinputformat.split.maxsize=10485760");
	arr_args.push_back("-Dmapred.input.fileinputformat.split.maxsize=10485760");
	arr_args.push_back("-Dmapreduce.input.fileinputformat.split.minsize=10485760");
	arr_args.push_back("-Dmapred.input.fileinputformat.split.minsize=10485760");*/

	/*arr_args.push_back("-D");
	arr_args.push_back("mapreduce.task.timeout=36000000");
	arr_args.push_back("-D");
	arr_args.push_back("mapred.max.split.size=10485760");
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.max.split.size=10485760");
	arr_args.push_back("-D");
	arr_args.push_back("mapred.min.split.size=0");
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.min.split.size=0");
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.input.fileinputformat.split.maxsize=10485760");
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.input.fileinputformat.split.maxsize=10485760");
	
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.input.fileinputformat.split.minsize=0");
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.input.fileinputformat.split.minsize=0");	

	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.job.maps=10");*/

	arr_args.push_back("-libjars");
	arr_args.push_back(fr_vars.hadoopgis_prefix + CUSTOM_JAR_REL_PATH);
	

	arr_args.push_back("-input");	
	arr_args.push_back(input_path);
	
	/*for(vector<string>::iterator it = input_paths.begin()+1; it != input_paths.end(); ++it) {
		arr_args.push_back("-input");
		arr_args.push_back(*it);
	}*/

	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + SKELETON);


	// should also work: 
	arr_args.push_back("-mapper");
	arr_args.push_back(SKELETON);
        //arr_args.push_back("cat");
	
	/*arr_args.push_back("-reducer");
	arr_args.push_back(SPACE_EXTRACTOR);

	// 2 reducers = 1 for outputting mbb and 1 for outputting spacial_dimension
	arr_args.push_back("-numReduceTasks");
	//arr_args.push_back("2");
	arr_args.push_back(fr_vars.numreducers_str);*/

	//arr_args.push_back("-numMapTasks");
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

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back("0");
	
	/*arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapred.max.split.size=10485760");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.max.split.size=10485760");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapred.min.split.size=0");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.min.split.size=0");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.input.fileinputformat.split.maxsize=10485760");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.input.fileinputformat.split.maxsize=10485760");
	
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.input.fileinputformat.split.minsize=0");
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.input.fileinputformat.split.minsize=0");*/

	cerr << "a: " << a  << " b: " << b << " c: " << c <<endl;
	cerr  <<"Optimal: " << optimal << endl;
	arr_args.push_back("-jobconf");
	stringstream ss;
	ss << "mapreduce.job.maps=" << optimal;
	arr_args.push_back(ss.str());
	
	//arr_args.push_back("-jobconf");
	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

#ifdef DEBUG
	cerr << "Extract 3D Skeleton program params: " << endl;
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in extracting Skeletons: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in extracting Skeletons: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}


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
	ss << COMPRESSION << " -o 2 " << fr_vars.sharedparams;
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
	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

#ifdef DEBUG
	cerr << "Compress data program params: " << endl;
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
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

bool build_voronoi(char *input, char *output, struct framework_vars &fr_vars) {

	vector<string> arr_args;
	arr_args.push_back(VORONOI);
	arr_args.push_back(input);
	arr_args.push_back(output);

	string programpath = fr_vars.hadoopgis_prefix + VORONOI;

#ifdef DEBUG
	cerr << programpath << endl;
	cerr << "Build Voronoi program params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in building Voronois: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in building Voronois: " << status << endl;
#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}



bool partition_data(string programpath, string input_path, 
	string output_path, string partitionmethod, int bucket_size, 
	string sharedparams, int step, double samplerate, struct framework_vars &fr_vars, 
	char *cachefilename, char *cachefilefullpath) {

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
		ss << MANIPULATE << " -o 1 " << sharedparams 
			<< " -c " << cachefilename ;
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
	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

#ifdef DEBUG
	cerr << "Partitioning params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
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
}


bool collect_stat(string programpath, string input_path, string output_path, 
	char *cachefilename, char *cachefilefullpath, struct framework_vars &fr_vars) {
	hdfs_delete(programpath, output_path);
	stringstream ss;
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};

	arr_args.push_back("-conf");
	arr_args.push_back("/home/vhoang/mapred-site.xml");
	
	arr_args.push_back("-input");
	arr_args.push_back(input_path);
	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + STAT_COLLECT_MAPPER);
	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + STAT_COLLECT_REDUCER);

	string strtmp(cachefilefullpath);
	arr_args.push_back("-file");
	arr_args.push_back(strtmp);

	arr_args.push_back("-mapper");
	ss.str("");
	ss << STAT_COLLECT_MAPPER << " -o 1 " << fr_vars.sharedparams 
			<< " -c " << cachefilename;
	arr_args.push_back(ss.str());

	arr_args.push_back("-reducer");
	arr_args.push_back(STAT_COLLECT_REDUCER);

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back(fr_vars.numreducers_str);

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");
	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

	#ifdef DEBUG
	cerr << "Collecting stats";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
	#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
		if (wait(&status)) {
			#ifdef DEBUG
			cerr << "Succeeded in collecting stats MBB: " << status << endl;
			#endif
		} else {
			#ifdef DEBUG
			cerr << "Failed in collecting stats " << status << endl;
			#endif
			exit(1);
		}
		return status == 0 ? true : false;
	}
}

bool duplicate_removal(string programpath, string input_path, string output_path, 
	struct framework_vars &fr_vars) {
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	arr_args.push_back("-input");
	arr_args.push_back(input_path);

	arr_args.push_back("-output");
	arr_args.push_back(output_path);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.hadoopgis_prefix + DUPLICATE_REMOVER);

	arr_args.push_back("-mapper");
	arr_args.push_back(DUPLICATE_REMOVER + " cat");

	//arr_args.push_back("-reducer None");
	arr_args.push_back("-reducer");
	arr_args.push_back(DUPLICATE_REMOVER + " uniq");
	// should also work: 
	// arr_args.push_back("cat");

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back(fr_vars.numreducers_str);
	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");
	arr_args.push_back("-cmdenv");
	arr_args.push_back(fr_vars.hadoopldlibpath);

#ifdef DEBUG
	cerr << "Removing duplicate program params: ";
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
#endif

	int status = 0;
	pid_t childpid;
	if (childpid = execute_command(programpath, arr_args)) {
		if (wait(&status)) {
#ifdef DEBUG
			cerr << "Succeeded in boundary handling: " << status << endl;
#endif
		} else {
#ifdef DEBUG
			cerr << "Failed in boundary handling: " << status << endl;
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


