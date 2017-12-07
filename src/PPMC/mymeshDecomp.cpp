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

#include "mymesh.h"
#include "frenetRotation.h"
#include "configuration.h"


/**
  * Start the next decompression operation.
  */
void MyMesh::startNextDecompresssionOp()
{
    if ((float)i_curOperationId / (i_nbQuantizations + i_nbDecimations) * 100 >= i_decompPercentage)
    {
        //std::cout << "End of mesh decompression." << std::endl;
        //std::cout << "Nb of LODs decompressed: " << i_curOperationId + 1 << std::endl;
        //std::cout << "First level to display: " << i_levelNotConvexId << std::endl;

        writeMeshOff(std::string(filePathOutput + ".off").c_str());

        for (MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
            hit->resetState();

        for (MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
            fit->resetState();

        operation = Idle;
        b_jobCompleted = true;
    }
    else
    {
        // Start the decoder.
        start_decoding(&rangeCoder);

        // Read the operation type.
        unsigned char i_operationType = decode_culshift(&rangeCoder, 1);
        decode_update(&rangeCoder, 1, i_operationType, 1 << 1);

        switch (i_operationType)
        {
        case DECIMATION_OPERATION_ID:
            beginUndecimationConquest();
            break;
        case QUANTIZATION_OPERATION_ID:
            beginAdaptiveUnquantization();
            break;
        default:
            assert(0);
            break;
        }
    }
}


/**
  * Begin an undecimation conquest.
  */
void MyMesh::beginUndecimationConquest()
{
    //printf("Begin undecimation conquest n°%u.\n", i_curDecimationId);

    for (MyMesh::Halfedge_iterator hit = halfedges_begin(); hit!=halfedges_end(); ++hit)
        hit->resetState();

    for (MyMesh::Face_iterator fit = facets_begin(); fit!=facets_end(); ++fit)
        fit->resetState();

    // Add the first halfedge to the queue.
    pushHehInit();

    operation = RemovedVertexCoding;
  //  printf("Removed vertex decoding begining.\n");

    f_avgSurfaceFaceWithCenterRemoved = 0;
    f_avgSurfaceFaceWithoutCenterRemoved = 0;
    i_nbFacesWithCenterRemoved = 0;
    i_nbFacesWithoutCenterRemoved = 0;
    i_nbGoodPredictions = 0;

    // Read if the prediction was used or not.
    unsigned test = decode_culshift(&rangeCoder, 1);
    b_predictionUsed = test;
    decode_update(&rangeCoder, 1, b_predictionUsed, 1 << 1);

    // Read the min values and the ranges.
    uint16_t i16_min;
    i16_min = decode_short(&rangeCoder);
    alphaBetaMin = *(int16_t *)&i16_min;
    unsigned alphaBetaRange = decode_short(&rangeCoder);

#ifdef USE_BIJECTION
    i16_min = decode_short(&rangeCoder);
    gammaMin = *(int16_t *)&i16_min;
    unsigned gammaRange = decode_short(&rangeCoder);
#endif

    // Init the range coder models.
    initqsmodel(&alphaBetaModel, alphaBetaRange, 18, 1 << 17, NULL, 0);
#ifdef USE_BIJECTION
    initqsmodel(&gammaModel, gammaRange, 18, 1 << 17, NULL, 0);
#endif
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 0);

    // Set the current operation.
    operation = UndecimationConquest;
}


/**
  * One undecimation step.
  */
void MyMesh::undecimationStep()
{
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();

        // If the face is already processed, pick the next halfedge:
        if (f->isConquered())
            continue;

        // Decode the face symbol.
        int syfreq, ltfreq;
        ltfreq = decode_culshift(&rangeCoder, 10);
        unsigned sym = qsgetsym(&connectModel, ltfreq);
        qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
        decode_update(&rangeCoder, syfreq, ltfreq, 1 << 10);
        // Update the model.
        qsupdate(&connectModel, sym);

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h;
        do
        {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isConquered())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        }
        while (hIt != h);

        bool b_split;
        float f_faceSurface = faceSurface(h);

        if (b_predictionUsed)
        {
            // Connectivity prediction.
            float f_nbConqueredFaces = i_nbFacesWithCenterRemoved + i_nbFacesWithoutCenterRemoved;
            if (b_useTriangleMeshConnectivityPredictionFaces && f_nbConqueredFaces != 0
                && fabs(i_nbFacesWithCenterRemoved / f_nbConqueredFaces
                        - i_nbFacesWithoutCenterRemoved / f_nbConqueredFaces) < 0.1)
            {
                // Triangle mesh connectivity prediction in the case the average area are similar.
                // Balance number between the splittable and non splittable neighbor faces.
                int i_neighborBalance = 0;

                // Compute the balance.
                Halfedge_around_facet_const_circulator hit = h->facet_begin(), hit_end = hit;
                CGAL_For_all(hit, hit_end)
                {
                    Face_const_handle fh = hit->opposite()->facet();

                    // If the current face is not already conquered,
                    // then it can't be taken into account for the prediction.
                    if (fh->isConquered())
                    {
                        if (fh->isSplittable())
                            i_neighborBalance++;
                        else
                            i_neighborBalance--;
                    }
                }

                // Determine if the face is split or not.
                if (sym == 1)
                    b_split = i_neighborBalance >= 0 ? false : true;
                else
                    b_split = i_neighborBalance >= 0 ? true : false;
            }
            else if ((f_avgSurfaceFaceWithCenterRemoved != 0
                 && fabs(f_faceSurface - f_avgSurfaceFaceWithCenterRemoved)
                    <= fabs(f_faceSurface - f_avgSurfaceFaceWithoutCenterRemoved))
                || f_avgSurfaceFaceWithoutCenterRemoved == 0)
            {
                // If the surface of the face is closer than the average surface
                // of the faces with a center vertex removed.
                b_split = sym == 1 ? true : false;
            }
            else
            {
                // If the surface of the face is closer than the average surface
                // of the faces without a center vertex removed.
                b_split = sym == 1 ? false : true;
            }
        }
        else
        {
            // No connectivity prediction.
            b_split = sym == 1 ? true : false;
        }

        // Update the average surfaces.
        updateAvgSurfaces(b_split, f_faceSurface);

        // Decode the geometry symbol.
        if (b_split)
            decodeGeometrySym(h, f);
        else
            f->setUnsplittable();

        return;
    }

   // printf("Removed vertex decoding completed.\n");

    // Stop the decoder.
    done_decoding(&rangeCoder);

    // Delete the models.
    deleteqsmodel(&connectModel);
    deleteqsmodel(&alphaBetaModel);
#ifdef USE_BIJECTION
    deleteqsmodel(&gammaModel);
#endif

    operation = InsertedEdgeDecoding;

    beginInsertedEdgeDecoding();
}


/**
  * Begin the inserted edge decoding conquest.
  */
void MyMesh::beginInsertedEdgeDecoding()
{
    // Add the first halfedge to the queue.
    pushHehInit();

    operation = InsertedEdgeDecoding;
  //  printf("Inserted edge decoding begining.\n");

    f_avgInsertedEdgesLength = 0;
    f_avgOriginalEdgesLength = 0;
    i_nbInsertedEdges = 0;
    i_nbOriginalEdges = 0;

    // Init the range coder models.
    initqsmodel(&connectModel, 2, 10, 1 << 9, NULL, 0);

    // Start the decoder.
    start_decoding(&rangeCoder);

    // Read if the prediction was used or not.
    b_predictionUsed = decode_culshift(&rangeCoder, 1);
    decode_update(&rangeCoder, 1, b_predictionUsed, 1 << 1);
}


/**
  * One step of the inserted edge coding conquest.
  */
void MyMesh::InsertedEdgeDecodingStep()
{
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        // Test if the edge has already been conquered.
        if (h->isProcessed())
            continue;

        // Mark the halfedge as processed.
        h->setProcessed();
        h->opposite()->setProcessed();

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h->next();
        while (hIt->opposite() != h)
        {
            if (!hIt->isProcessed() && !hIt->isNew())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }

        assert(!hIt->isNew());

        bool b_original = true;
        float f_edgeLen = edgeLen(h);

        // Test if there is a symbol for this edge.
        // There is no symbol if the two faces of an egde are unsplitable.
        if (h->facet()->isSplittable()
            || h->opposite()->facet()->isSplittable())
        {
            // Decode the edge symbol.
            int syfreq, ltfreq;
            ltfreq = decode_culshift(&rangeCoder, 10);
            unsigned sym = qsgetsym(&connectModel, ltfreq);
            qsgetfreq(&connectModel, sym, &syfreq, &ltfreq);
            decode_update(&rangeCoder, syfreq, ltfreq, 1 << 10);
            // Update the model.
            qsupdate(&connectModel, sym);

            // Determine if the edge is original or not.

            if (b_predictionUsed)
            {
                if ((f_avgOriginalEdgesLength != 0
                     && fabs(f_edgeLen - f_avgOriginalEdgesLength)
                     <= fabs(f_edgeLen - f_avgInsertedEdgesLength))
                        || f_avgInsertedEdgesLength == 0)
                {
                    b_original = sym == 1 ? true : false;
                }
                else
                {
                    b_original = sym == 1 ? false : true;
                }
            }
            else
                b_original = sym == 0 ? true : false;

            // Mark the edge to be removed.
            if (!b_original)
                h->setAdded();
        }

        // Update the average edge lengths.
        updateAvgEdgeLen(b_original, f_edgeLen);

        return;
    }

    //printf("Inserted edge decoding completed.\n");

    // Stop the decoder.
    done_decoding(&rangeCoder);

    // Delete the models.
    deleteqsmodel(&connectModel);

    if (b_useLiftingScheme)
        lift(true); // Unlift the vertex positions.

    insertRemovedVertices();

    removeInsertedEdges();

    //printf("Number of vertices: %lu - Number of faces: %lu\n",
    //       size_of_vertices(), size_of_facets());

    if (i_mode == DECOMPRESSION_MODE_WRITE_ALL_ID)
        writeCurrentOperationMesh(filePathOutput, i_curOperationId + 1);

    i_curDecimationId++; // Increment the current decimation operation id.
    i_curOperationId++;
    operation = Idle;
}


/**
  * Insert center vertices.
  */
void MyMesh ::insertRemovedVertices()
{
   // printf("Insert removed vertices.\n");

    // Add the first halfedge to the queue.
    pushHehInit();

    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        gateQueue.pop();

        Face_handle f = h->facet();

        // If the face is already processed, pick the next halfedge:
        if (f->isProcessed())
            continue;

        // Mark the face as processed.
        f->setProcessedFlag();

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h;
        do
        {
            Halfedge_handle hOpp = hIt->opposite();
            assert(!hOpp->is_border());
            if (!hOpp->facet()->isProcessed())
                gateQueue.push(hOpp);
            hIt = hIt->next();
        }
        while (hIt != h);

        assert(!h->isNew());

        if (f->isSplittable())
        {
            // Insert the vertex.
            Point p = getPos(getQuantizedPos(barycenter(h)) + f->getResidual());
            Halfedge_handle hehNewVertex = create_center_vertex(h);
            hehNewVertex->vertex()->point() = p;

            // Mark all the created edges as new.
            Halfedge_around_vertex_circulator Hvc = hehNewVertex->vertex_begin();
            Halfedge_around_vertex_circulator Hvc_end = Hvc;
            CGAL_For_all(Hvc, Hvc_end)
            {
                Hvc->setNew();
                Hvc->opposite()->setNew();
            }
        }
    }
}


/**
  * Remove all the marked edges.
  */
void MyMesh::removeInsertedEdges()
{
    //printf("Remove inserted edges.\n");

    for (MyMesh::Halfedge_iterator hit = halfedges_begin();
         hit!=halfedges_end(); ++hit)
    {
        if(hit->isAdded())
            join_facet(hit);
    }
}


/**
  * Decode the geometry symbols.
  */
void MyMesh::decodeGeometrySym(Halfedge_handle heh_gate, Face_handle fh)
{
#ifdef USE_BIJECTION
    Vector t1 = CGAL::NULL_VECTOR;
    Vector t2 = CGAL::NULL_VECTOR;

    Vector normal = computeNormal(heh_gate);
    if (normal == CGAL::NULL_VECTOR)
    {
        t1 = Vector(1,0,0);
        t2 = Vector(0,1,0);
        normal = Vector(0,0,1);
    }
    else
        determineFrenetFrame(heh_gate, normal, t1, t2);
#endif

    int coord[3];
    for (unsigned i = 0; i < 3; ++i)
    {
        int syfreq, ltfreq;
#ifdef USE_BIJECTION
        if (i < 2)
        {
#endif
            // Decode the alpha and beta symbols.
            ltfreq = decode_culshift(&rangeCoder, 18);
            unsigned sym = qsgetsym(&alphaBetaModel, ltfreq);
            qsgetfreq(&alphaBetaModel, sym, &syfreq, &ltfreq);
            decode_update(&rangeCoder, syfreq, ltfreq, 1 << 18);
            // Update the alpha and beta model.
            qsupdate(&alphaBetaModel, sym);
            // Store the value.
            coord[i] = alphaBetaMin + sym;
#ifdef USE_BIJECTION
        }
        else
        {
            // Encode the gamma symbol.
            ltfreq = decode_culshift(&rangeCoder, 18);
            unsigned sym = qsgetsym(&gammaModel, ltfreq);
            qsgetfreq(&gammaModel, sym, &syfreq, &ltfreq);
            decode_update(&rangeCoder, syfreq, ltfreq, 1 << 18);
            // Update the gamma model.
            qsupdate(&gammaModel, sym);
            // Store the value.
            coord[i] = gammaMin + sym;
        }
#endif
    }

#ifdef USE_BIJECTION
    VectorInt frenetCoord(coord[0], coord[1], coord[2]);
    // Compute the barycenter with the same way as decoding to avoid computation errors.
    VectorInt correction = invFrenetRotation(frenetCoord, t1, t2, normal);
#else
    VectorInt correction(coord[0], coord[1], coord[2]);
    // Compute the barycenter with the same way as decoding to avoid computation errors.
#endif

    if (b_useCurvaturePrediction)
        correction = correction - avgLaplacianVect(heh_gate) / INV_ALPHA;

    fh->setSplittable();
    fh->setResidual(correction);
}


/**
  * Begin the adaptive unquantization operation.
  */
void MyMesh::beginAdaptiveUnquantization()
{
    //printf("Adaptive unquantization begining.\n");

    for(MyMesh::Vertex_iterator vit = vertices_begin(); vit!=vertices_end(); ++vit)
          vit->resetState();

    // Add the first halfedge to the queue.
    pushHehInit();

    operation = AdaptiveUnquantization;

    // Init the range coder models.
    initqsmodel(&quantModel, 8, 12, 2000, NULL, 0);
}


/**
  * One step of the inserted edge coding conquest.
  */
void MyMesh::adaptiveUnquantizationStep()
{
    while (!gateQueue.empty())
    {
        Halfedge_handle h = gateQueue.front();
        Vertex_handle vh = h->vertex();
        gateQueue.pop();

        // Test if the edge has already been conquered.
        if (vh->isConquered())
            continue;

        // Mark the vertex as processed.
        vh->setConquered();

        // Add the other halfedges to the queue
        Halfedge_handle hIt = h->next();
        while (hIt->opposite() != h)
        {
            if (!hIt->vertex()->isConquered())
                gateQueue.push(hIt);
            hIt = hIt->opposite()->next();
        }

        // Decode the vertex symbol.
        int syfreq, ltfreq;
        ltfreq = decode_culshift(&rangeCoder, 12);
        unsigned sym = qsgetsym(&quantModel, ltfreq);
        qsgetfreq(&quantModel, sym, &syfreq, &ltfreq);
        decode_update(&rangeCoder, syfreq, ltfreq, 1 << 12);
        // Update the model.
        qsupdate(&quantModel, sym);

        std::map<unsigned, unsigned> cellMap = determineCellSymbols(h, false);
        unsigned i_cellId = cellMap[sym];

        VectorInt correction;
        switch (i_cellId)
        {
        case 0:
            correction = VectorInt(0, 0, 0);
            break;
        case 1:
            correction = VectorInt(1, 0, 0);
            break;
        case 2:
            correction = VectorInt(0, 1, 0);
            break;
        case 3:
            correction = VectorInt(1, 1, 0);
            break;
        case 4:
            correction = VectorInt(0, 0, 1);
            break;
        case 5:
            correction = VectorInt(1, 0, 1);
            break;
        case 6:
            correction = VectorInt(0, 1, 1);
            break;
        case 7:
            correction = VectorInt(1, 1, 1);
            break;
        }

        PointInt oldPosQuant = getQuantizedPos(vh->point());
        PointInt p = oldPosQuant + VectorInt(CGAL::ORIGIN, oldPosQuant) + correction; // = 2 * oldPosQuant + correction

        vh->point() = Point((p.x() + 0.5) * f_quantStep * (1 << i_curQuantizationId - 1) + bbMin.x(),
                            (p.y() + 0.5) * f_quantStep * (1 << i_curQuantizationId - 1) + bbMin.y(),
                            (p.z() + 0.5) * f_quantStep * (1 << i_curQuantizationId - 1) + bbMin.z());

        return;
    }

    // Stop the decoder.
    done_decoding(&rangeCoder);

    // Delete the model.
    deleteqsmodel(&quantModel);

    if (i_mode == DECOMPRESSION_MODE_WRITE_ALL_ID)
        writeCurrentOperationMesh(filePathOutput, i_curOperationId + 1);

    i_curQuantizationId--;
    i_curOperationId++;

    //printf("Adaptive unquantization completed.\n");

    operation = Idle;
}
