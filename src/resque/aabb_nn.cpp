/*
 * aabb_nn.cpp
 *
 *  Created on: Sep 24, 2019
 *      Author: teng
 */

#include <resque/resque_3d.hpp>
using namespace std;
using namespace SpatialIndex;

/*
 * perform nearest neighbor query on data sets with AABB tree
 * the kernel for CGAL we used in this module is Simple Cartesian
 * */
int join_bucket_nn_rtree(struct query_op &stop, struct query_temp &sttemp) {

	assert(stop.join_cardinality == 2 && "cannot conduct nn on same data set");

	/* Indicates where original data is mapped to */
	int idx1 = SID_1;
	int idx2 = SID_2;
	int pairs = 0; // number of satisfied results
	double low[3], high[3];  // Temporary value placeholders for MBB
	int kneighs = 2; // kNN

	int len1 = sttemp.mbbdata[idx1].size();
	int len2 = sttemp.mbbdata[idx2].size();
 	if (len1 <= 0 || len2 <= 0) {
		 return 0;
	}

	try {

		/* Handling for special nearest neighbor query */
		// build the actual spatial index for input polygons from idx2
		IStorageManager *storage = NULL;
		ISpatialIndex *spidx = NULL;
		if (! build_index_geoms(sttemp.mbbdata[idx2], spidx, storage)) {
			#ifdef DEBUG
			cerr << "Building index on geometries from set 2 has failed" << endl;
			#endif
			return -1;
		}
#ifdef DEBUG
		else{
			cerr << "index is built on data set 2"<<endl;
		}
#endif

		//vector<vector<long>> id1_offset2; // for each nuclei id, the vector of its nearest blood vessels' offset
		unordered_set<int> unique_nn_id2; // the unique set of nearest blood vessels' offset
		unordered_map<int, vector<int>> nn_id2; // mapping between the unique offset and length
		vector<Sc_Point> nuclei_pts;

		vector<struct mbb_3d *> geom_mbb1 = sttemp.mbbdata[idx1];
		for (int i = 0; i < geom_mbb1.size(); i++) {

			struct mbb_3d * env1 = geom_mbb1[i];

			low[0] = env1->low[0];
			low[1] = env1->low[1];
			low[2] = env1->low[2];
			high[0] = env1->high[0];
			high[1] = env1->high[1];
			high[2] = env1->high[2];

			// Temporary value placeholders for MBB
			double np[3];
			MyVisitor vis;
			/* R-tree intersection check */

			np[0] = (low[0]+high[0])*0.5;
			np[1] = (low[1]+high[1])*0.5;
			np[2] = (low[2]+high[2])*0.5;
			SpatialIndex::Point nuclei_centroid(SpatialIndex::Point(np, 3));
			vis.matches.clear();
			/* Find kNN objects*/
			spidx->nearestNeighborQuery(1, nuclei_centroid, vis);
			#ifdef DEBUG
			cerr << "mbb data size " << sttemp.mbbdata[idx2].size() << TAB << " and found " << vis.matches.size() << endl;
			#endif
			for (uint32_t j = 0; j < vis.matches.size(); j++) {
				#ifdef DEBUG
				cerr << "Query the nearest neighbor between " << i << TAB << vis.matches[j] << endl;
				#endif

				// push the offset and length of its nearest blood vessels
				//long offset = sttemp.offsetdata[idx2][vis.matches[j]], length = sttemp.lengthdata[idx2][vis.matches[j]];
				nn_id2[i].push_back(vis.matches[j]);
				// record the unique blood vessels' info
				unique_nn_id2.insert(vis.matches[j]);
			}
			nuclei_pts.push_back(Sc_Point(np[0], np[1], np[2]));
		}

		#ifdef DEBUG
		cerr << "Unique NN poly is: " << unique_nn_id2.size() << endl;
		#endif


		/* for each unique nearest blood vessel, construct the AABB tree*/
		unordered_map<int, Sc_Tree*> id2_aabbtree; // map between unique id of blood vessel and its AABB tree
		// for each mentioned object in data set 2, build an AABB tree
		for(auto it = unique_nn_id2.begin(); it != unique_nn_id2.end(); ++it ){

			long offset = sttemp.offsetdata[idx2][*it];
			long length = sttemp.lengthdata[idx2][*it];

			Sc_Polyhedron geom2 = sc_extract_geometry(offset, length, stop.decomp_lod, stop, sttemp, 1);
//			Sc_Tree *facet_tree = new Sc_Tree(); // zero-constructed
//			int valid = 0;
//			int invalid = 0;
//			for(Sc_Polyhedron::Facet_iterator fit = geom2.facets_begin(),
//					end = geom2.facets_end(); fit!=end; ++fit) {
//				Sc_Polyhedron::Point a(fit->halfedge()->vertex()->point()),
//									 b(fit->halfedge()->next()->vertex()->point()),
//									 c(fit->halfedge()->prev()->vertex()->point());
//
//				if(!CGAL::collinear(a,b,c)) {
//					//facet_tree->insert(Facet_primitive(fit, *poly_item->polyhedron(),
//					facet_tree->insert(Sc_Primitive(fit, geom2, geom2));
//					valid++;
//				} else {
//					invalid++;
//				}
//			}
//			cerr << "Valid and invalid :" << valid << TAB << invalid << endl;
//      		facet_tree->build();
//			facet_tree->accelerate_distance_queries();
//			id2_aabbtree[*it] = facet_tree;
			try{
				Sc_Tree *tree = new Sc_Tree(faces(geom2).first, faces(geom2).second, geom2);
				id2_aabbtree[*it] = tree;
			}catch(std::exception &exc){
				cerr << exc.what()<<endl;
			}
		}
		#ifdef DEBUG
		cerr << "Done creating AABB tree" << endl;
		#endif

		/* for each nuclei, calculate distance by searching the AABB tree of its k nearest blood vessels*/
		Sc_Point nuclei_point;
		Sc_Tree *aabbtree;  // AABB tree for each blood vessel
		int pairs = 0; // number of satisfied results
		for (int j = 0; j < nuclei_pts.size(); j++) {
			//cerr << vis.matches[j] << endl;
			//poly = poly2_geom[poly_index[j]];
			//cerr << *poly << endl;
			nuclei_point = nuclei_pts[j];
			vector<int> ids = nn_id2[j];
			#ifdef DEBUG
			cerr << "# of NN vessel is: " << ids.size() << endl;
			#endif
			//TODO which one should we use?
			//for(int m = 0; m < 2; m++){
			for(int m = 0; m < ids.size(); m++){
				aabbtree = id2_aabbtree[ids[m]];
				aabbtree->accelerate_distance_queries();
				#ifdef DEBUG
				cerr << "Checking distance calc between " << j << TAB << ids[m] << endl;
				cerr << "point  coords " << nuclei_point << endl;
				#endif
				//double squared_dist = CGAL::to_double(CGAL::squared_distance(nnp, nuclei_point));
				cerr<<"getting distance"<<endl;

				Sc_FT sqd = aabbtree->squared_distance(nuclei_point);
				cerr<<"got distance"<<endl;
				sttemp.nn_distance = sqrt((double)CGAL::to_double(sqd));


				cout <<  j << TAB << nuclei_point.x() << TAB << nuclei_point.y() << TAB << nuclei_point.z() << TAB << sttemp.nn_distance << endl;
				pairs++;
			}
		}

		cerr << "Done with tile" << endl;

		delete spidx;
		delete storage;

		// release aabb tree of blood vessels
		for (auto it = id2_aabbtree.begin(); it != id2_aabbtree.end(); ++it ) {
			delete it->second;
		}
		id2_aabbtree.clear();

	} catch (Tools::Exception& e) {
		std::cerr << "******ERROR******" << std::endl;
		#ifdef DEBUG
		cerr << e.what() << std::endl;
		#endif
		return -1;
	} // end of catch

	return pairs ;
}
