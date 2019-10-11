
#include <progparams/queryprocessor_params_3d.hpp>

namespace po = boost::program_options;
using namespace std;

inline void remove_slash(string &str){
	if (str.at(str.size() - 1) == '/') {
		str = str.substr(0, str.size() - 1);
	}
}

bool extract_params(int argc, char **argv, struct framework_vars &fr_vars) {
	string binpath;
	string querytype;
	try {
		po::options_description desc("Options");
		desc.add_options()
			// common parameters
			("help,h", "This help message")
			("querytype,q", po::value<string>(&querytype), "Query type [ partition | compress | join]")
			("binpath",po::value<string>(&binpath), "path to the binary executables")
			("outputpath,o", po::value<string>(&fr_vars.output_path), "Output path")
			("overwrite", "Overwrite existing hdfs directories")
			("numreducers,n", po::value<int>(&fr_vars.numreducers), "The number of reducers")

			("input1,a", po::value<string>(&fr_vars.input_path_1), "HDFS file path to data set 1")
			("input2,b", po::value<string>(&fr_vars.input_path_2), "HDFS file path to data set 2")
			// partition
			("samplingrate,s", po::value<double>(&fr_vars.sampling_rate), "Sampling rate (0, 1]")
			("partitioner,t", po::value<string>(&fr_vars.partition_method), "Partitioning method ([fg_3d | ot_3d ]")
			("bucket,k", po::value<long>(&fr_vars.bucket_size), "Fine-grain level tile size for spjoin")

			// resque
			("lod,l", po::value<int>(&fr_vars.decomp_lod) , "Decompression LOD. (0, 100]. Default is 100.")
			("distance,d", po::value<double>(&fr_vars.distance), "Distance (used for certain predicates)")
			("predicate,p", po::value<string>(&fr_vars.predicate), "Predicate for spatial join and nn queries "
					"[ st_intersects | st_touches | st_crosses | st_contains | st_adjacent | st_disjoint "
					"| st_equals | st_dwithin | st_within | st_overlaps | st_nn_voronoi | st_nn_rtree ] ")
			("compressed_data_path",po::value<string>(&fr_vars.compressed_data_path), "path to the combined compressed spatial data")
			;
		po::variables_map vm;        
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);   

		if (vm.count("help")) {
			cerr << desc << endl;
			return false;
		}
		if (!vm.count("querytype")) {
			cerr << desc << endl;
			cerr << "operation type must be given" << endl;
			return false;
		}
		fr_vars.query_type = get_operation(querytype);
		if(!vm.count("binpath")){
			cerr << desc << endl;
			cerr << "path to the binary executables must be given"<<endl;
			return false;
		}
		if(!vm.count("outputpath")){
			cerr << desc << endl;
			cerr << "path to the output folder must be given"<<endl;
			return false;
		}

		remove_slash(fr_vars.output_path);
		if (vm.count("overwrite")) {
			fr_vars.overwritepath = true;
		}

		// some general constrains
		switch(fr_vars.query_type){
		case JOIN:
			if(!vm.count("predicate") || get_join_predicate(fr_vars.predicate.c_str())==ST_ERROR){
				cerr << desc << endl;
				cerr << "type of predict should be set properly"<<endl;
				return false;
			}
			break;
		case COMPRESS:
			if(!vm.count("input1")||!vm.count("input2")){
				cerr << desc << endl;
				cerr << "path for input1 and input2 should be given"<<endl;
				return false;
			}
			remove_slash(fr_vars.input_path_1);
			remove_slash(fr_vars.input_path_2);
			break;
		case PARTITION:
			// sample rate, bucket size, partitioner can be specified
			if (fr_vars.partition_method != PARTITION_FG_3D
				&& fr_vars.partition_method != PARTITION_OT_3D) {
				cerr << desc << endl;
				cerr << "Invalid partitioner." << endl;
				return false;
			}
			break;
		case DUPLICATE_REMOVAL:
			break;
		default:
			break;
		}

		/* Update environment variables with HADOOP_HOME defined*/
		fr_vars.hadoopcmdpath = getHadoopCmdPath();
		fr_vars.hdfscmdpath = getHdfsCmdPath();
		fr_vars.streaming_path = getHadoopJarPath();

		fr_vars.binary_prefix = binpath + SLASH;
		#ifdef DEBUG
		cerr << "Hadoop commands path:" << fr_vars.hadoopcmdpath  << endl;
		cerr << "Number reducers: " << fr_vars.numreducers << endl;
		#endif

	} catch (exception& e) {
		cerr << "error here: " << e.what() << "\n";
		return false;
	} catch (...) {
		cerr << "Exception of unknown type!\n";
		return false;
	}

	std::stringstream tmpss;

	// updating the paths for output folders
	tmpss.str("");
	tmpss << fr_vars.output_path << "_mbb";
	fr_vars.mbb_output = tmpss.str();

	tmpss.str("");
	tmpss << fr_vars.output_path << "_inputresque";
	fr_vars.resque_input = tmpss.str();

	tmpss.str("");
	tmpss << fr_vars.output_path << "_partidx";;
	fr_vars.partitionpath = tmpss.str();

	tmpss.str("");
	tmpss << fr_vars.partitionpath << "/part*";
	fr_vars.partitionpathout = tmpss.str();
	
	tmpss.str("");
	tmpss << fr_vars.output_path << "_joinout";
	fr_vars.joinoutputpath = tmpss.str();

	return true;
}

