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

//#include "PPMC/mymesh.h"
//#include "PPMC/configuration.h"
//#include "PPMC/frenetRotation.h"
//#include <unistd.h>
//#include <algorithm>


MyMesh::MyMesh(char filename[], 
	       //std::string filePathOutput,
               unsigned i_decompPercentage,
               const int i_mode,
               unsigned i_quantBits,
               bool b_useAdaptiveQuantization,
               bool b_useLiftingScheme,
               bool b_useCurvaturePrediction,
               bool b_useConnectivityPredictionFaces,
               bool b_useConnectivityPredictionEdges,
               bool b_allowConcaveFaces,
               bool b_useTriangleMeshConnectivityPredictionFaces,
	       std::string & ss,
		char* offset,
		long length,
		char *buffer_loc
) :
    CGAL::Polyhedron_3< CGAL::Simple_cartesian<float>, MyItems >(), i_mode(i_mode),
    b_jobCompleted(false), operation(Idle),
    i_curDecimationId(0), i_curQuantizationId(0), i_curOperationId(0),
    i_levelNotConvexId(0), b_testConvexity(false), connectivitySize(0),
    geometrySize(0), i_quantBits(i_quantBits), dataOffset(0),
   // filePathOutput(filePathOutput),
    i_decompPercentage(i_decompPercentage),
    osDebug(&fbDebug),
    b_useAdaptiveQuantization(b_useAdaptiveQuantization),
    b_useLiftingScheme(b_useLiftingScheme),
    b_useCurvaturePrediction(b_useCurvaturePrediction),
    b_useConnectivityPredictionFaces(b_useConnectivityPredictionFaces),
    b_useConnectivityPredictionEdges(b_useConnectivityPredictionEdges),
    b_useTriangleMeshConnectivityPredictionFaces(b_useTriangleMeshConnectivityPredictionFaces),
    ss(ss), offset(offset), length(length)
{
    // Create the compressed data buffer.
    //p_data = new char[BUFFER_SIZE];
    p_data = buffer_loc;
    // Fill the buffer with 0.
 
    for (size_t i = 0; i < BUFFER_SIZE; ++i) {
       p_data[i] = 0;
   }

    // Initialise the range coder structure.
    rangeCoder.p_data = p_data;
    rangeCoder.p_dataOffset = &dataOffset;

    //std::cout << "p_data: " << (void *)p_data << std::endl;
    //std::cout << "dataOffset: " << dataOffset << std::endl;
    if (i_mode == COMPRESSION_MODE_ID) // Compression mode.
    {
	
        //std::filebuf fb;
        //fb.open (filename, std::ios::in);
        //if(fb.is_open())
        //{
            //std::istream is(&fb);
	    std::istringstream is;
	    is.str(ss);
            //std::cerr << ss << std::endl;
            is >> *this;
         //   std::cerr << is.size() << std::endl;

	    
            //fb.close();

            if (keep_largest_connected_components(1) != 0)
            {
                std::cout << "Can't compress the mesh." << std::endl;
                std::cout << "The codec doesn't handle meshes with several connected components." << std::endl;
                exit(EXIT_FAILURE);
            }

            if (!is_closed())
            {
                std::cout << "Can't compress the mesh." << std::endl;
                std::cout << "The codec doesn't handle meshes with borders." << std::endl;
                exit(EXIT_FAILURE);
            }

            /* The special connectivity prediction scheme for triangle
               mesh is not used if the current mesh is not a pure triangle mesh. */
            if (!is_pure_triangle())
                this->b_useTriangleMeshConnectivityPredictionFaces = false;

            computeBoundingBox();
            determineQuantStep();
            quantizeVertexPositions();

#if 0
            // Output the initial quantified mesh in an off file.
            writeMeshOff("mesh_quant.off");
#endif

      /*      printf("Bounding box min coordinates: %f %f %f.\n",
                   bbMin.x(), bbMin.y(), bbMin.z());
            printf("Quantization step: %f\n", f_quantStep); */

            // Set the vertices of the edge that is the departure of the coding and decoding conquests.
            vh_departureConquest[0] = halfedges_begin()->opposite()->vertex();
            vh_departureConquest[1] = halfedges_begin()->vertex();

            if (!b_allowConcaveFaces) {
                b_testConvexity = true;
	    }
        //}
    }
    else // Decompression mode.
    {
        //readCompressedFile(filename);
	//std::cerr << "reading decompressed shm" << std::endl;
	
	//std::cerr << "dataOffset" << dataOffset << "\t" << (long) p_data  << std::endl;
	readCompressedFile(offset, length);
        // std::cerr << "dataOffset" << dataOffset << "\t" << (long) p_data  << std::endl;
	//std::cerr << "reading base mesh" << std::endl;
        readCompressedData();
	//std::cerr << "done reading bash mesh" << std::endl;
        // std::cerr << "dataOffset" << dataOffset << "\t" << (long) p_data  << std::endl;


	//sleep(5);
	//std::cerr << "done reading bash mesh" << std::endl;
        if (i_mode == DECOMPRESSION_MODE_WRITE_ALL_ID){
	   // std::cerr << "i_mode == DECOMPRESSION_MODE_WRITE_ALL_ID" << std::endl;
            writeCurrentOperationMesh(filePathOutput, 0);
	}
	//else
	//    std::cerr << "i_mode != DECOMPRESSION_MODE_WRITE_ALL_ID" << std::endl;
        
	// Set the vertices of the edge that is the departure of the coding and decoding conquests.
        vh_departureConquest[0] = vertices_begin();
        vh_departureConquest[1] = ++vertices_begin();
	//std::cerr << "vertices_begin(): " << vh_departureConquest[0] << std::endl;
	//std::cerr << "++vertices_begin(): " << ++vertices_begin() << std::endl;
	
	//std::cerr << "before batch operation" << std::endl;
        if (i_levelNotConvexId - i_nbDecimations > 0)
        {
            // Decompress until the first convex LOD is reached.
     //       std::cerr << "happening here" << std::endl;
            while (i_curDecimationId < i_levelNotConvexId) {
       //     	std::cerr << "first batch" << std::endl;
                batchOperation();
            }
	    
         //   std::cerr << "second batch" << std::endl;
            batchOperation();
        }
    }
    //std::cerr << "done batch operation" << std::endl;
    i_nbVerticesInit = size_of_vertices();
    i_nbFacetsInit = size_of_facets();
    //printf("Number of vertices: %lu.\nNumber of faces: %lu.\n",
    //      i_nbVerticesInit, i_nbFacetsInit);

    //fbDebug.open("debug.txt", std::ios::out | std::ios::trunc);
}


MyMesh::~MyMesh()
{
    //delete[] p_data;
   // fbDebug.close();
}



/**
  * Perform one step of the current operation.
  */
void MyMesh::stepOperation()
{
    if (b_jobCompleted)
        return;
    switch (operation)
    {
    case Idle:
        if (i_mode == COMPRESSION_MODE_ID)
            startNextCompresssionOp();
        else
            startNextDecompresssionOp();
        break;
    case DecimationConquest:
        decimationStep();
        break;
    case RemovedVertexCoding:
        RemovedVertexCodingStep();
        break;
    case InsertedEdgeCoding:
        InsertedEdgeCodingStep();
        break;
    case AdaptiveQuantization:
        adaptiveQuantizationStep();
        break;
    case UndecimationConquest:
        undecimationStep();
        break;
    case InsertedEdgeDecoding:
        InsertedEdgeDecodingStep();
        break;
    case AdaptiveUnquantization:
        adaptiveUnquantizationStep();
        break;
    default:
        break;
    }
}


/**
  * Perform one batch of steps of the current operation.
  */
void MyMesh::batchOperation()
{
    if (b_jobCompleted)
        return;
    switch (operation)
    {
    case Idle:
        if (i_mode == COMPRESSION_MODE_ID)
            startNextCompresssionOp();
        else
            startNextDecompresssionOp();
        break;
    case DecimationConquest:
        while(operation == DecimationConquest)
            decimationStep();
        break;
    case RemovedVertexCoding:
        while(operation == RemovedVertexCoding)
            RemovedVertexCodingStep();
        break;
    case InsertedEdgeCoding:
        while(operation == InsertedEdgeCoding)
            InsertedEdgeCodingStep();
        break;
    case AdaptiveQuantization:
        while(operation == AdaptiveQuantization)
            adaptiveQuantizationStep();
        break;
    case UndecimationConquest:
        while(operation == UndecimationConquest)
            undecimationStep();
        break;
    case InsertedEdgeDecoding:
        while(operation == InsertedEdgeDecoding)
            InsertedEdgeDecodingStep();
        break;
    case AdaptiveUnquantization:
        while (operation == AdaptiveUnquantization)
            adaptiveUnquantizationStep();
        break;
    default:
        break;
    }
}


/**
  * Finish completely the current operation.
  */
void MyMesh::completeOperation()
{
    while (!b_jobCompleted)
        batchOperation();
}


// Compute the mesh vertex bounding box.
void MyMesh::computeBoundingBox()
{
    //printf("Compute the mesh bounding box.\n");
    std::list<Point> vList;
    for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
        vList.push_back(vit->point());
    MyKernel::Iso_cuboid_3 bBox = CGAL::bounding_box(vList.begin(), vList.end());

    bbMin = bBox.min();
    bbMax = bBox.max();

    bbMin0 = bbMin;
    bbMax0 = bbMax;

    Vector bbVect = bbMax - bbMin;
    f_bbVolume = bbVect.x() * bbVect.y() * bbVect.z();
}


// Get the bounding box diagonal length.
float MyMesh::getBBoxDiagonal() const
{
    return sqrt(Vector(bbMin, bbMax).squared_length());
}


/**
  * Get the bounding box center point
  */
Vector MyMesh::getBBoxCenter() const
{
    return ((bbMax - CGAL::ORIGIN)
            + (bbMin - CGAL::ORIGIN)) / 2;
}


/**
  * Determine the quantization step.
  */
void MyMesh::determineQuantStep()
{
    //printf("Determine the quantization step.\n");

    float f_maxRange = 0;
    for (unsigned i = 0; i < 3; ++i)
    {
        float range = bbMax[i] - bbMin[i];
        if (range > f_maxRange)
            f_maxRange = range;
    }
    f_quantStep = f_maxRange / (1 << i_quantBits);

    if (b_useLiftingScheme)
        bbMin = bbMin - Vector(1, 1, 1) * (1 << (i_quantBits - 1)) * f_quantStep;

    if (b_useAdaptiveQuantization)
        f_adaptQuantRescaling = 10 / f_maxRange;
}


// Compute and store the quantized positions of the mesh vertices.
void MyMesh::quantizeVertexPositions()
{
    //printf("Quantize the vertex positions.\n");
    unsigned i_maxCoord = 1 << i_quantBits;

    // Update the positions to fit the quantization.
    for (MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
    {
        Point p = vit->point();
        PointInt quantPoint = getQuantizedPos(p);

        if (!b_useLiftingScheme)
        {
            // Make sure the last value is in the range.
            assert(quantPoint.x() <= i_maxCoord);
            assert(quantPoint.y() <= i_maxCoord);
            assert(quantPoint.z() <= i_maxCoord);
            /* The max value is the unique that have to to be reassigned
               because else it would get out the range. */
            quantPoint = PointInt(quantPoint.x() == i_maxCoord ? i_maxCoord - 1 : quantPoint.x(),
                                  quantPoint.y() == i_maxCoord ? i_maxCoord - 1 : quantPoint.y(),
                                  quantPoint.z() == i_maxCoord ? i_maxCoord - 1 : quantPoint.z());
        }

        Point newPos = getPos(quantPoint);
        vit->point() = newPos;
    }
}


/**
  * Quantize a position
  */
PointInt MyMesh::getQuantizedPos(Point p) const
{
    return PointInt((p.x() - bbMin.x()) / (f_quantStep * (1 << i_curQuantizationId)),
                    (p.y() - bbMin.y()) / (f_quantStep * (1 << i_curQuantizationId)),
                    (p.z() - bbMin.z()) / (f_quantStep * (1 << i_curQuantizationId)));
}


/**
  * Get a position from the quantized coordinates.
  */
Point MyMesh::getPos(PointInt p) const
{
    return Point((p.x() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.x(),
                 (p.y() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.y(),
                 (p.z() + 0.5) * f_quantStep * (1 << i_curQuantizationId) + bbMin.z());
}

void MyMesh::printPdata()
{
     for (size_t myidx = 0; myidx < 1000; ++myidx) {
       //printf("%02x ", p_data[myidx]);
       printf("%x ", (unsigned) p_data[myidx]);

     }
     printf("\n\n");

}
