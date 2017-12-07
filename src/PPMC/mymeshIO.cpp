/*****************************************************************************
* Copyright (C) 2011 Adrien Maglo
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

//#include <fstream>

//#include "mymesh.h"
//#include "mymeshBaseBuilder.h"
//#include "configuration.h"

// #include <unistd.h>
/**
  * Write the compressed data to the buffer.
  */
void MyMesh::writeCompressedData()
{
    i_nbDecimations = i_curDecimationId + 1;
    i_nbQuantizations = i_curQuantizationId;

    // Write the base mesh.
    writeBaseMesh();
    //printf("Base mesh size: %lu bytes.\n", (size_t)ceil((connectivitySize + geometrySize) / 8.0));

    // Write the number of bytes needed to decode the base mesh.
    osDebug << connectivitySize + geometrySize
            /*<< "," << connectivitySize << "," << geometrySize*/ << std::endl;

    //printf("Writing the compressed data.\n");

    unsigned i_deci = i_curDecimationId;
    unsigned i_quant = i_curQuantizationId - 1;

 //   std::cerr << "Params: " << i_deci << "\t" << i_quant << std::endl;

    while (!typeOfOperation.empty())
    {
        unsigned i_curOperationType = typeOfOperation.back();
        typeOfOperation.pop_back();

        switch (i_curOperationType)
        {
        case DECIMATION_OPERATION_ID:
            encodeRemovedVertices(i_deci);
            encodeInsertedEdges(i_deci);
            i_deci--;
            break;
        case QUANTIZATION_OPERATION_ID:
            encodeAdaptiveQuantization(adaptiveQuantSym[i_quant]);
            i_quant--;
            break;
        default:
            assert(0);
            break;
        }



        // Write the number of bytes needed to decode the current LOD.
        osDebug << connectivitySize + geometrySize /*<< ","
                << connectivitySize << "," << geometrySize*/ << std::endl;
    }

    //printf("Connectivity size: %lu bytes -> %.2f bits per vertex.\n",
    //       (size_t)ceil(connectivitySize / 8.0), connectivitySize / (float)i_nbVerticesInit);
    //printf("Geometry size: %lu bytes -> %.2f bits per vertex.\n",
    //       (size_t)ceil(geometrySize / 8.0), geometrySize / (float)i_nbVerticesInit);

    size_t totalSize = connectivitySize + geometrySize;
    //printf("Total size: %lu bytes -> %.2f bits per vertex.\n",
    //       (size_t)ceil(totalSize / 8.0), totalSize / (float)i_nbVerticesInit);

    //printf("Id of the first not convex decimation: %d.\n", i_levelNotConvexId);

    //std::cout << "totalSize: " << totalSize << std::endl;
}


/**
  * Read the compressed data from the buffer.
  * \return the bounding box min position.
  */
void MyMesh::readCompressedData()
{
    //printf("Got here: readCompressedData()");
    // Read the base mesh.
    readBaseMesh();

 //   printf("Bounding box min coordinates: %f %f %f.\n",
 //         bbMin.x(), bbMin.y(), bbMin.z());
 //    printf("Quantization step: %f\n", f_quantStep);
}


// Write a given number of bits in a buffer.
void writeBits(uint32_t data, unsigned i_nbBits, char *p_dest,
               unsigned &i_bitOffset, size_t &offset)
{
    assert(i_nbBits <= 25);

    uint32_t dataToAdd = data << (32 - i_nbBits - i_bitOffset);
    // Swap the integer bytes because the x86 architecture is little endian.
    dataToAdd = __builtin_bswap32(dataToAdd); // Call a GCC builtin function.

    // Write the data.
    *(uint32_t *)p_dest |= dataToAdd;

    // Update the size and offset.
    offset += (i_bitOffset + i_nbBits) / 8;
    i_bitOffset = (i_bitOffset + i_nbBits) % 8;
}


/**
  * Read a given number of bits in a buffer.
  */
uint32_t readBits(unsigned i_nbBits, char *p_src,
                  unsigned &i_bitOffset, size_t &offset)
{
    assert(i_nbBits <= 25);

    // Build the mask.
    uint32_t mask = 0;
    for (unsigned i = 0; i < 32 - i_bitOffset; ++i)
        mask |= 1 << i;
    // Swap the mask bytes because the x86 architecture is little endian.
    mask = __builtin_bswap32(mask); // Call a GCC builtin function.

    uint32_t data = *(uint32_t *)p_src & mask;

    // Swap the integer bytes because the x86 architecture is little endian.
    data = __builtin_bswap32(data); // Call a GCC builtin function.

    data >>= 32 - i_nbBits - i_bitOffset;

    // Update the size and offset.
    offset += (i_bitOffset + i_nbBits) / 8;
    i_bitOffset = (i_bitOffset + i_nbBits) % 8;

    return data;
}


// Write a floating point number in the data buffer.
void MyMesh::writeFloat(float f)
{
    *(float *)(p_data + dataOffset) = f;
    dataOffset += sizeof(float);
}


/**
  * Read a floating point number in the data buffer.
  */
float MyMesh::readFloat()
{
    float f = *(float *)(p_data + dataOffset);
    dataOffset += sizeof(float);
    return f;
}


// Write a 16 bits integer in the data buffer
void MyMesh::writeInt16(int16_t i)
{
    *(int16_t *)(p_data + dataOffset) = i;
    dataOffset += sizeof(int16_t);
}


/**
  * Read a 16 bits integer in the data buffer.
  */
int16_t MyMesh::readInt16()
{
    int16_t i = *(int16_t *)(p_data + dataOffset);
    dataOffset += sizeof(int16_t);
    return i;
}


// Write the base mesh.
void MyMesh::writeBaseMesh()
{
   // printf("Writing the base mesh.\n");
   // printf("before writing basemesh");
   // printPdata();

    // Write the bounding box min coordinate.
    for (unsigned i = 0; i < 3; ++i)
        writeFloat(bbMin[i]);
    // Write the quantization step.
    writeFloat(f_quantStep);

    geometrySize += 4 * sizeof(float) * 8;

    unsigned i_nbVerticesBaseMesh = size_of_vertices();
    unsigned i_nbFacesBaseMesh = size_of_facets();
    unsigned i_nbBitsPerVertex = ceil(log(i_nbVerticesBaseMesh) / log(2));

    unsigned i_bitOffset = 0;
    dataOffset++;
 
    //printf("after quantization step");
    //printPdata();
 

    // Write the codec option status.
    writeBits(b_useAdaptiveQuantization, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(b_useLiftingScheme, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(b_useCurvaturePrediction, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(b_useConnectivityPredictionFaces, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(b_useConnectivityPredictionEdges, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(b_useTriangleMeshConnectivityPredictionFaces, 1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    geometrySize += 3;
    connectivitySize += 4;
   
   // printf("After codec writing ");
   // printPdata();

    // Write the geometry quantization of the mesh.
    assert(i_quantBits - 1 < 1 << 4);
    writeBits(i_quantBits - 1, 4, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    geometrySize += 4;

    // Write the number of level of decimations.
    assert(i_nbDecimations < 1 << 6);
    writeBits(i_nbDecimations, 6, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    connectivitySize += 6;

    // Write the number of adaptive quantizations.
    assert(i_nbQuantizations < 1 << 6);
    writeBits(i_nbQuantizations, 6, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    connectivitySize += 6;

    // Write the number of non-convex level of details.
    writeBits(i_nbDecimations - i_levelNotConvexId, 6, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    connectivitySize += 6;
    //printf("After writing number of LoDs");
    //printPdata();

    // Write the number of vertices and faces on 16 bits.
    //printf("Base mesh: %u vertices and %u faces.\n", i_nbVerticesBaseMesh, i_nbFacesBaseMesh);
    writeBits(i_nbVerticesBaseMesh, 16, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    writeBits(i_nbFacesBaseMesh, 16, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    connectivitySize += 32;
   // printf("After writing # of vertices and faces");
   // printPdata();

    // Write the base mesh vertex coordinates.
    unsigned i_nbAdditionalBitsGeometry = b_useLiftingScheme ? LIFTING_NB_ADDITIONAL_BITS_GEOMETRY : 0;

    // Write the vertices of the edge that is the departure of the coding conquests.
    for (unsigned j = 0; j < 2; ++j)
    {
        PointInt p = getQuantizedPos(vh_departureConquest[j]->point());
        for (unsigned i = 0; i < 3; ++i)
        {
            assert(p[i] < 1 << i_quantBits - i_curQuantizationId + i_nbAdditionalBitsGeometry);
            writeBits(p[i], i_quantBits - i_curQuantizationId + i_nbAdditionalBitsGeometry,
                      p_data + dataOffset - 1, i_bitOffset, dataOffset);
        }
        vh_departureConquest[j]->setId(j);
    }

    // Write the other vertices.
    size_t id = 2;
    for (MyMesh::Vertex_iterator vit = vertices_begin();
         vit != vertices_end(); ++vit)
    {
        if (vit == vh_departureConquest[0] || vit == vh_departureConquest[1])
            continue;

        PointInt p = getQuantizedPos(vit->point());

        // Write the coordinates.
        for (unsigned i = 0; i < 3; ++i)
        {
            assert(p[i] < 1 << i_quantBits - i_curQuantizationId + i_nbAdditionalBitsGeometry);
            writeBits(p[i], i_quantBits - i_curQuantizationId + i_nbAdditionalBitsGeometry, p_data + dataOffset - 1, i_bitOffset, dataOffset);
        }
        // Set an id to the vertex.
        vit->setId(id++);
    }
    geometrySize += i_nbVerticesBaseMesh * 3 * (i_quantBits - i_curQuantizationId);

    // Write the base mesh face vertex indices.
    for (MyMesh::Facet_iterator fit = facets_begin();
         fit != facets_end(); ++fit)
    {
        unsigned i_faceDegree = fit->facet_degree();
        unsigned i_code = i_faceDegree - 3;
        assert(i_code < 1 << NB_BITS_FACE_DEGREE_BASE_MESH);
        writeBits(i_code, NB_BITS_FACE_DEGREE_BASE_MESH, p_data + dataOffset - 1, i_bitOffset, dataOffset);

        Halfedge_around_facet_const_circulator hit(fit->facet_begin()), end(hit);
        do
        {
            // Write the current vertex id.
            writeBits(hit->vertex()->getId(), i_nbBitsPerVertex, p_data + dataOffset - 1, i_bitOffset, dataOffset);
        }
        while(++hit != end);

        connectivitySize += i_nbBitsPerVertex * i_faceDegree + NB_BITS_FACE_DEGREE_BASE_MESH;
    }

    if(i_bitOffset == 0)
        dataOffset--;

//     printf("Final pdata: ");
//     printPdata();

}


// Read the base mesh.
void MyMesh::readBaseMesh()
{
    //printf("Reading the base mesh.\n");

    // Read the bounding box min coordinate.
    float coord[3];
    for (unsigned i = 0; i < 3; ++i)
        coord[i] = readFloat();
    bbMin = Point(coord[0], coord[1], coord[2]);

    // Read the quantization step.
    f_quantStep = readFloat();

    unsigned i_bitOffset = 0;
    dataOffset++;

    // Read the codec option status.
    b_useAdaptiveQuantization = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    b_useLiftingScheme = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    b_useCurvaturePrediction = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    b_useConnectivityPredictionFaces = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    b_useConnectivityPredictionEdges = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    b_useTriangleMeshConnectivityPredictionFaces = readBits(1, p_data + dataOffset - 1, i_bitOffset, dataOffset);

    // Read the geometry quantization of the mesh.
    i_quantBits = readBits(4, p_data + dataOffset - 1, i_bitOffset, dataOffset) + 1;

    // Read the number of level of detail.
    i_nbDecimations = readBits(6, p_data + dataOffset - 1, i_bitOffset, dataOffset);

    // Read the number of quantization operations.
    i_nbQuantizations = i_curQuantizationId = readBits(6, p_data + dataOffset - 1, i_bitOffset, dataOffset);

    // Read the number of non convex level of details.
    i_levelNotConvexId = readBits(6, p_data + dataOffset - 1, i_bitOffset, dataOffset);

    // Set the mesh bounding box.
    unsigned i_nbQuantStep = 1 << i_quantBits;
    bbMax = bbMin + Vector(i_nbQuantStep * f_quantStep,
                           i_nbQuantStep * f_quantStep,
                           i_nbQuantStep * f_quantStep);

    unsigned i_nbVerticesBaseMesh = readBits(16, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    unsigned i_nbFacesBaseMesh = readBits(16, p_data + dataOffset - 1, i_bitOffset, dataOffset);
    unsigned i_nbBitsPerVertex = ceil(log(i_nbVerticesBaseMesh) / log(2));

    std::deque<Point> *p_pointDeque = new std::deque<Point>();
    std::deque<uint32_t *> *p_faceDeque = new std::deque<uint32_t *>();
    unsigned i_nbAdditionalBitsGeometry = b_useLiftingScheme ? LIFTING_NB_ADDITIONAL_BITS_GEOMETRY : 0;

    // Read the vertex positions.
    for (unsigned i = 0; i < i_nbVerticesBaseMesh; ++i)
    {
        uint32_t p[3];
        for (unsigned j = 0; j < 3; ++j)
            p[j] = readBits(i_quantBits - i_curQuantizationId + i_nbAdditionalBitsGeometry,
                            p_data + dataOffset - 1, i_bitOffset, dataOffset);
        PointInt posInt(p[0], p[1], p[2]);
        Point pos = getPos(posInt);
        p_pointDeque->push_back(pos);
    }

    // Read the face vertex indices.
    for (unsigned i = 0; i < i_nbFacesBaseMesh; ++i)
    {
        uint32_t *f = new uint32_t[(1 << NB_BITS_FACE_DEGREE_BASE_MESH) + 3];

        // Write in the first cell of the array the face degree.
        f[0] = readBits(NB_BITS_FACE_DEGREE_BASE_MESH, p_data + dataOffset - 1, i_bitOffset, dataOffset) + 3;

        for (unsigned j = 1; j < f[0] + 1; ++j)
            f[j] = readBits(i_nbBitsPerVertex, p_data + dataOffset - 1, i_bitOffset, dataOffset);

        p_faceDeque->push_back(f);
    }

    // Let the builder do its job.
    MyMeshBaseBuilder<HalfedgeDS> builder(p_pointDeque, p_faceDeque);
    delegate(builder);

    // Free the memory.
    for (unsigned i = 0; i < p_faceDeque->size(); ++i)
        delete[] p_faceDeque->at(i);
    delete p_faceDeque;
    delete p_pointDeque;

    if(i_bitOffset == 0)
        dataOffset--;

    //printf("Done with Reading the base mesh in the function.\n");
}


// Write the compressed data from the buffer to a file.
int MyMesh::writeCompressedFile() const
{
    int i_ret = 1;

    //std::cout << "Write the compressed file " << filePathOutput << "." << std::endl;

    std::filebuf fb;
    fb.open(filePathOutput.c_str(), std::ios::out | std::ios::trunc);
     //std::cout << "p_data: " << (void *)p_data << std::endl;
    if (fb.is_open())
    {
        if (fb.sputn(p_data, dataOffset) == (std::streamsize)dataOffset)
            i_ret = 0;
        fb.close();
	//std::cout <<"dataOffset: " << dataOffset << std::endl;
	//std::cout << "p_data: " << (void *)p_data << std::endl;
    }

    return i_ret;
}

/*
// Read the compressed data from the file and store them in a buffer.
int MyMesh::readCompressedFile(char psz_filePath[])
{
    int i_ret = 1;

    //printf("Read the compressed file '%s'.\n", psz_filePath);

    std::filebuf fb;
    fb.open(psz_filePath, std::ios::in);

    if (fb.is_open())
    {
        std::streamsize dataSize = fb.in_avail();
        if (fb.sgetn(p_data, dataSize) == (std::streamsize)dataSize)
            i_ret = 0;
        fb.close();
    }

    return i_ret;
}
*/

/*
// Write the compressed data from the buffer to a file.
int MyMesh::writeCompressedFile() const
{
    int i_ret = 1;

    std::cout << "Write the compressed file " << filePathOutput << "." << std::endl;

    std::filebuf fb;
    fb.open(filePathOutput.c_str(), std::ios::out | std::ios::trunc);
     std::cout << "p_data: " << (void *)p_data << std::endl;
    if (fb.is_open())
    {
        if (fb.sputn(p_data, dataOffset) == (std::streamsize)dataOffset)
            i_ret = 0;
        fb.close();
	std::cout <<"dataOffset: " << dataOffset << std::endl;
	std::cout << "p_data: " << (void *)p_data << std::endl;
    }

    return i_ret;
}
*/

// Read the compressed data from the file and store them in a buffer.
int MyMesh::readCompressedFile(char* offset, long length)
{
    /*int i_ret = 1;

    printf("Read the compressed file '%s'.\n", psz_filePath);

    std::filebuf fb;
    fb.open(psz_filePath, std::ios::in);

    if (fb.is_open())
    {
        std::streamsize dataSize = fb.in_avail();
        if (fb.sgetn(p_data, dataSize) == (std::streamsize)dataSize)
            i_ret = 0;
        fb.close();
    }

    return i_ret;*/

/*
  //std::cerr << "got here and doing memcpy " <<  std::endl;
  std::cerr << "got here and doing memcpy " <<  ( (long) offset) << "\t" << length << "\t" <<  ( (long ) p_data)  <<  std::endl;
  
    
    for (size_t i = 0; i < 1000; ++i) {
	printf("%02X ", offset[i]);
    }
   printf("\n");
   for (size_t i = 0; i < 1000; ++i) {
	printf("%02X ", p_data[i]);
	}
   printf("\n");
*/


    // Copy data from shm to p_data
    memcpy(p_data, offset, length);


  /*
    for (size_t i = 0; i < 1000; ++i){
	printf("%02X ", p_data[i]);
    }
   printf("\n");
*/
 //
 // Hoang commented out this
 //dataOffset = 0;
 //
 //

    //std::cout <<"dataOffset: " << dataOffset << std::endl;
    //	std::cout << "p_data: " << (void *)p_data << std::endl;
    //dataSize = length;
    return 0;
}


// Write the mesh in an off file.
void MyMesh::writeMeshOff(const char psz_filePath[]) const
{
    std::filebuf fb;
    fb.open(psz_filePath, std::ios::out | std::ios::trunc);
    if(fb.is_open())
    {
        std::ostream os(&fb);
        os << *this;
    }
}


void MyMesh::writeCurrentOperationMesh(std::string pathPrefix, unsigned i_id) const
{
    // Output the current mesh in an off file.
    std::ostringstream fileName;
    fileName << pathPrefix << "_" << i_id << ".off";
    writeMeshOff(fileName.str().c_str());
}
