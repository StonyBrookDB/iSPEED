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
	if (!extract_params(argc, argv, fr_vars)) {
		return 1;
	}
	
#ifdef DEBUG
	time_t start_exec_time, end_exec_time;
	time(&start_exec_time);
#endif

	/*
	 *
	 * process operations
	 *
	 * */

	switch(fr_vars.query_type){
	case COMPRESS:
		execute_compress(fr_vars);
		break;
	case PARTITION:
		execute_partition(fr_vars);
		break;
	case JOIN:
		execute_spjoin(fr_vars);
		break;
	case DUPLICATE_REMOVAL:
		execute_duplicate_removal(fr_vars);
		break;
	default:
		break;
	}

#ifdef DEBUG
	time(&end_exec_time);
	double total_exec_time = difftime(end_exec_time,start_exec_time);
	cerr << "********************************************" << endl;
	cerr << "Total execution time: " << total_exec_time << " seconds." << endl;
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

bool compress_data(struct framework_vars &fr_vars) {
	hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.mbb_output);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	
	arr_args.push_back("-D");
	arr_args.push_back("mapreduce.task.timeout=36000000");
	
	arr_args.push_back("-input");
	arr_args.push_back(fr_vars.input_path_1);
	if (fr_vars.join_cardinality > 1) {
		arr_args.push_back("-input");
		arr_args.push_back(fr_vars.input_path_2);
	}

	arr_args.push_back("-output");
	arr_args.push_back(fr_vars.mbb_output);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.binary_prefix + COMPRESSION);

	//ppmc is used to compress each object and extract mbbs
	//the path for the second dataset is used for identifying
	//which dataset each mapper processed
	arr_args.push_back("-mapper");
	stringstream ss;
	ss << COMPRESSION << " " << fr_vars.input_path_2;
	arr_args.push_back(ss.str());
	ss.str("");
	arr_args.push_back("-numReduceTasks");
	arr_args.push_back("0");

	#ifdef DEBUG
	cerr << "Compress data program params: " << endl;
	for(vector<string>::iterator it = arr_args.begin(); it != arr_args.end(); ++it) {
		cerr << *it << " ";
	}
	cerr << endl;
	#endif

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(fr_vars.hadoopcmdpath, arr_args))) {
		if (wait(&status)) {
			cerr << "Succeeded in compressing data: " << status << endl;
		} else {
			cerr << "Failed in compressing data: " << status << endl;
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}


bool partition_data(struct framework_vars &fr_vars) {

	hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.partitionpath);
	stringstream ss;
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};

	arr_args.push_back("-input");
	arr_args.push_back(fr_vars.resque_input);
	arr_args.push_back("-output");
	arr_args.push_back(fr_vars.partitionpath);

	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.binary_prefix + MBB_SAMPLER);
	arr_args.push_back("-file");
	arr_args.push_back(fr_vars.binary_prefix + fr_vars.partition_method);

	arr_args.push_back("-mapper");
	ss.str("");
	ss << MBB_SAMPLER << " " << fr_vars.sampling_rate;
	arr_args.push_back(ss.str());

	arr_args.push_back("-reducer");
	ss.str("");
	long bucket_size = max(static_cast<int>(floor(fr_vars.bucket_size * fr_vars.sampling_rate)), 1);
	ss  << fr_vars.partition_method << " -b " << bucket_size
		<<" --min_x "<<fr_vars.spinfo.space_low[0]
		<<" --min_y "<<fr_vars.spinfo.space_low[1]
		<<" --min_z "<<fr_vars.spinfo.space_low[2]
		<<" --max_x "<<fr_vars.spinfo.space_high[0]
		<<" --max_y "<<fr_vars.spinfo.space_high[1]
		<<" --max_z "<<fr_vars.spinfo.space_high[2];

	arr_args.push_back(ss.str());

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back("1");
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
	if ((childpid = execute_command(fr_vars.hadoopcmdpath, arr_args))) {
		if (wait(&status)) {
			cerr << "Succeeded in partitioning " << status << endl;
		} else {
			cerr << "Failed in partitioning: " << status << endl;
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}

bool join_data(struct framework_vars &fr_vars, char *cachefilefullpath) {
	char *cachefilename = strrchr(cachefilefullpath, '/'); // pointing to just the name of the cache file
	cachefilename++; // Advance the pointer past the slash delimiter character

	hdfs_delete(fr_vars.hadoopcmdpath, fr_vars.joinoutputpath);
	vector<string> arr_args = {"hadoop", "jar", fr_vars.streaming_path};
	arr_args.push_back("-files");
	stringstream ss;
	ss << fr_vars.binary_prefix + MANIPULATE << ","
	   << fr_vars.binary_prefix + RESQUE << ","
	   << cachefilefullpath;
	arr_args.push_back(ss.str());

	arr_args.push_back("-input");
	arr_args.push_back(fr_vars.resque_input);

	arr_args.push_back("-output");
	arr_args.push_back(fr_vars.joinoutputpath);

	// the mapper phase assign each object to
	// different tiles with its mbb
	arr_args.push_back("-mapper");
	ss.str("");
	ss << MANIPULATE << " " << cachefilename;
	arr_args.push_back(ss.str());

	// the resque tool received mbbs for spatial join
	arr_args.push_back("-reducer");
	ss.str("");
	ss <<RESQUE
	   <<" -l "<<fr_vars.decomp_lod
	   <<" -j "<<fr_vars.join_cardinality
	   <<" -p "<<fr_vars.predicate
	   <<" -s "<<fr_vars.size_of_compressed_data;
	arr_args.push_back(ss.str()); // Offset to account for tile id and join index

	arr_args.push_back("-numReduceTasks");
	arr_args.push_back(to_string(fr_vars.numreducers));

	arr_args.push_back("-jobconf");
	arr_args.push_back("mapreduce.task.timeout=36000000");

	int status = 0;
	pid_t childpid;
	if ((childpid = execute_command(fr_vars.hadoopcmdpath, arr_args))) {
		if (wait(&status)) {
			cerr << "Succeeded in sp join: " << status << endl;
		} else {
			cerr << "Failed in sp join: " << status << endl;
			exit(1);
		}
		return status == 0 ? true : false;
	}
	return false;
}


inline void profile_raw_data(struct framework_vars &fr_vars){
	// Find the size of datasets 1 and 2
	fr_vars.size_1 = hdfs_get_size(fr_vars.hadoopcmdpath, fr_vars.input_path_1);
	fr_vars.spinfo.total_size = fr_vars.size_1;

	// If there is a 2nd dataset, also check if it has been loaded and obtain its size
	if (fr_vars.join_cardinality > 1) {
		// loaded_2 = hdfs_check_data(hadoopcmdpath, input_path_2 + "/" + PARTITION_FILE_NAME);
		fr_vars.size_2 = hdfs_get_size(fr_vars.hadoopcmdpath, fr_vars.input_path_2);
		if (fr_vars.size_2 >= 0) {
			fr_vars.spinfo.total_size += fr_vars.size_2;
		}
	}

	cerr << "Total size of 1 " << fr_vars.size_1 << endl;
	if (fr_vars.join_cardinality > 1) {
		cerr << "Total size of 2 " << fr_vars.size_2 << endl;
	}
	cerr << "Total object size: " << fr_vars.spinfo.total_size << endl;
}

inline void profile_compressed_data(struct framework_vars &fr_vars){
	// we assume the combiner already been executed and the
	// space information can be retrieved from the combined binary file.
	struct stat results;
	if (stat(fr_vars.compressed_data_path.c_str(), &results) != 0 ||
			results.st_size<(6*sizeof(double)+sizeof(long))){
		std::cerr<<"error reading compressed file "<<fr_vars.compressed_data_path<<std::endl;
		exit(0);
	}
	fr_vars.size_of_compressed_data = results.st_size;
	int fd = open(fr_vars.compressed_data_path.c_str(), O_RDONLY);
	lseek(fd,results.st_size-6*sizeof(double)-sizeof(long),SEEK_CUR);
	read(fd,fr_vars.spinfo.space_low,3*sizeof(double));
	read(fd,fr_vars.spinfo.space_high,3*sizeof(double));
	read(fd,&fr_vars.spinfo.num_objects,sizeof(long));
	close(fd);
}


/*
 * 1: in the first phase, load data from all data sets, and
 *    extract mbbs and compress them with ppmc
 * */
void execute_compress(struct framework_vars &fr_vars){

	cerr << "\nExecuting compressing\n"<< endl;
	profile_raw_data(fr_vars);

	// compress the input data
	// -- Extract object MBB and grab space dimension
	// -- and compress data at the same time to binary files
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.mbb_output)) {
		cerr << "\nCompressing & Extracting MBBs\n" << endl;
		if (!compress_data(fr_vars)) {
			cerr << "Failed extracting MBB"  << endl;
			exit(1);
		}
	}
	cerr<<"done compressing and extracting MBBs"<<endl;
}

/*
 * 2: in the second phase Run combiner on primary node and then loaders on all nodes
 *    to load the compressed data into a shared memory
 *
 * */


/*
 * 3: in the third phase, partition the space with the mbbs get in the
 * 	  first and second phase. partition index files are generated
 * 	  which will be used in the last phase
 *
 * */
void execute_partition(struct framework_vars &fr_vars){
	// profile the compressed data for more information
	// used for partitioning the space
	profile_compressed_data(fr_vars);

	//partition objects into tiles
	cerr << "\nExecuting partitioning\n" << endl;
	// Setting default block size if it has not been set
	if (fr_vars.bucket_size < 0) {
		// Bucket size was not set
		double blockSize = 16000000; // approximately 16MB
		fr_vars.bucket_size = max(static_cast<int>(floor(
				blockSize/fr_vars.size_of_compressed_data*fr_vars.spinfo.num_objects)), 1);
	}
	cerr << "Bucket size: " << fr_vars.bucket_size << endl;
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.partitionpath)) {
		if (!partition_data(fr_vars)) {
			cerr << "Failed partitioning" << endl;
			exit(1);
		}
	}
	cerr << "done partitioning space" <<endl;
}

/*
 * 4: do the real spatial join in the fourth phase
 *
 * */
void execute_spjoin(struct framework_vars &fr_vars) {
	profile_compressed_data(fr_vars);

	cerr << "\nExecuting spatial joins\n" << endl;
	//generate a temporary file for storing the partition information
	int tmpfd = mkstemp(nametemplate);
	char *tmpFile = nametemplate;
	close(tmpfd);
	tmpfd = open(tmpFile, O_RDWR | O_CREAT | O_TRUNC , 0777);
	if(!hdfs_cat(fr_vars.hadoopcmdpath, fr_vars.partitionpathout, tmpfd)){
		cerr<<"cannot load partition index from path "<<fr_vars.partitionpathout<<endl;
		exit(1);
	}
	close(tmpfd);
	cerr << "Temp file name to hold partition boundary: " << tmpFile << endl;

	// Most important-heavy-lifting work here
	// inputresque: is the list of all objects MBBs with offset and length of the level 0 compression
	//  e.g. 1203 231 242341 12424 12441 42414 0 200 (starts at 0 byte and 200 bytes)
	// Spatial join step
	if (fr_vars.overwritepath || !hdfs_check_data(fr_vars.hadoopcmdpath, fr_vars.joinoutputpath)) {
		if (!join_data(fr_vars,tmpFile)) {
			cerr << "Failed spatial join" << endl;
			exit(1);
		}
	}
#ifndef DEBUG
	remove(tmpFile);
#endif
	cerr << "Done with spatial join." << endl;

	cout.flush();
	cerr.flush();
}

void execute_duplicate_removal(struct framework_vars &fr_vars){
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
}




