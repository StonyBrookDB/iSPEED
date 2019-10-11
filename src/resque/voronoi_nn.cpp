/*
 * voronoi_nn.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: teng
 */

#include "resque_3d.hpp"
using namespace std;

/*
 * in each tile, each large object like vessel is extracted into
 * skeletons. Then a voronoi graph is built on top of those skeletons
 * and objects in data set 1 will get there distance to the nearest object
 * in data set 2 by looking up the voronoi graph
 *
 * */
int join_bucket_nn_voronoi(struct query_op &stop, struct query_temp &sttemp) {

	/* Indicates where original data is mapped to */
	int idx1 = SID_1;
	int idx2 = SID_2;
	int nuclei_id = 0;
	double low[3], high[3];  // Temporary value placeholders for MBB

	try {
		// extract the geometry from dataset2 (compressed blood vessels) and extract skeleton
		Sc_Skeleton *skeleton = NULL;

		std::vector<Sc_Point> P;
		vector<struct mbb_3d *> geom_mbb2 = sttemp.mbbdata[idx2];
#ifdef DEBUG
		std::cerr<<"size of data set 2: "<<geom_mbb2.size()<<std::endl;
#endif
		for (int i = 0; i < geom_mbb2.size(); i++) {
			//use the advanced way to extract skeleton, the simple one
			//has bugs on extracting skeleton from polyhedron compressed
			//by PPMC
		  	extract_skeleton_advance(sttemp.offsetdata[idx2][i],
		  			sttemp.lengthdata[idx2][i], stop.decomp_lod, P);
		}
		// building their Delaunay triangulation (Voronoi).
		Sc_Delaunay T(P.begin(), P.end());
#ifdef DEBUG
		std::cerr<<"voronoi graph is built"<<std::endl;
#endif

		// For each nuclei, find its nearest blood vessel by checking voronoi
		vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];

		for (int i = 0; i < geom_mbb1.size(); i++) {

			struct mbb_3d * env1 = geom_mbb1[i];
			Sc_Point nuclei_centroid((env1->low[0]+env1->high[0])*0.5,
					(env1->low[1]+env1->high[1])*0.5, (env1->low[2]+env1->high[2])*0.5);

			Sc_Point nnp = T.nearest_vertex(nuclei_centroid)->point();
			double squared_dist = CGAL::to_double(CGAL::squared_distance(nnp, nuclei_centroid));
			double distance = sqrt(squared_dist);
#ifdef DEBUG
			std::cerr<<"distance: "<<distance<<std::endl;
#endif
			cout << nuclei_id << TAB << nuclei_centroid.x() << TAB << nuclei_centroid.y()
					<< TAB << nuclei_centroid.z() << TAB << distance << "\n";
			nuclei_id++;
		}

	} catch (Tools::Exception& e) {
		std::cerr << "******ERROR******" << std::endl;
		#ifdef DEBUG
		cerr << e.what() << std::endl;
		#endif
		return -1;
	} // end of catch

	return nuclei_id ;
}

