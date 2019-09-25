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

	assert(stop.join_cardinality == 2 && "cannot do self nearest neighbor query");
	/* Indicates where original data is mapped to */
	int idx1 = SID_1;
	int idx2 = SID_2;

	int nuclei_id = 0;
	double low[3], high[3];  // Temporary value placeholders for MBB

	try {

		int len1 = sttemp.mbbdata[idx1].size();
		int len2 = sttemp.mbbdata[idx2].size();

	 	if (len1 <= 0 || len2 <= 0) {
			 return 0;
		}

		// extract the geometry from dataset2 (compressed blood vessels) and extract skeleton
		Sc_Skeleton *skeleton = NULL;

		std::vector<Sc_Point> P;
		vector<struct mbb_3d *> geom_mbb2 = sttemp.mbbdata[idx2];

		for (int i = 0; i < geom_mbb2.size(); i++) {

			#ifdef DEBUG
			// Extract geometry from compressed data in the shared memeory segment
			cerr << "obj: " << i << endl;
			cerr << "offset: " << sttemp.offsetdata[idx2][i] << endl;
			cerr << "length: " << sttemp.lengthdata[idx2][i] << endl;
			#endif

			try {
				// extract the skeleton of input polyhedron
#ifdef DEBUG
		  		cerr << "extracting the Skeleton!" << endl;
#endif
		  		Sc_Skeleton skeleton = extract_skeleton(sttemp.offsetdata[idx2][i],
		  				sttemp.lengthdata[idx2][i], stop.decomp_lod);
		  		std::cerr << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
		  		std::cerr << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
		  	//	// Output all the edges of the skeleton.
		  	//	BOOST_FOREACH(Sc_Skeleton_edge e, edges(skeleton)){
		  	//		const Point& s = skeleton[source(e, skeleton)].point;
		  	//		const Point& t = skeleton[target(e, skeleton)].point;
		  	//		std::cerr << "2 "<< s << " " << t << "\n";
		  	//	}
//		  		Sc_Skeleton skeleton;
//		  		Sc_Polyhedron geom = sc_extract_geometry(sttemp.offsetdata[idx2][i],
//		  				sttemp.lengthdata[idx2][i], stop.decomp_lod, stop, sttemp, 1);
//		  		CGAL::extract_mean_curvature_flow_skeleton(geom, skeleton);

		  		BOOST_FOREACH(Sc_Skeleton_vertex v, vertices(skeleton)){
					Sc_Point p = skeleton[v].point;
					P.push_back(p);
					cerr << p<< endl;
					//P.push_back(CGAL_Point3(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z())));
				}
#ifdef DEBUG
		  		cerr << "extracted one Skeleton!" << endl;
#endif
			} catch (const std::exception &exc) {
				cerr << "******Extract Skeleton Error******" << endl;
				cerr << exc.what() << endl;
				return -1;
			}
		}
		// building their Delaunay triangulation (Voronoi).
		Sc_Delaunay T(P.begin(), P.end());
#ifdef DEBUG
		std::cerr<<" voronoi graph is built"<<std::endl;
#endif

		// For each nuclei, find its nearest blood vessel by checking voronoi
		vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];

		for (int i = 0; i < geom_mbb1.size(); i++) {

			struct mbb_3d * env1 = geom_mbb1[i];
			Sc_Point nuclei_centroid((env1->low[0]+env1->high[0])*0.5,
					(env1->low[1]+env1->high[1])*0.5, (env1->low[2]+env1->high[2])*0.5);

			Sc_Point nnp = T.nearest_vertex(nuclei_centroid)->point();
			double squared_dist = CGAL::to_double(CGAL::squared_distance(nnp, nuclei_centroid));
			sttemp.nn_distance = sqrt(squared_dist);

			cout << nuclei_id << TAB << nuclei_centroid.x() << TAB << nuclei_centroid.y()
					<< TAB << nuclei_centroid.z() << TAB << sttemp.nn_distance << "\n";

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

