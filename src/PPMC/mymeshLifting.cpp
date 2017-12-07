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
#include "configuration.h"


/**
  * Lift the vertex positions.
  */
void MyMesh::lift(bool b_unlift)
{
    /*if (b_unlift)
        printf("Unlift.\n");
    else
        printf("Lift.\n");*/

    for (MyMesh::Vertex_iterator vit = vertices_begin();
         vit != vertices_end(); ++vit)
          vit->resetState();

    // Add the first halfedge to the queue.
    pushHehInit();

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

        //for each of the neighbour vertices v_i (including the removed vertices of the neighbouring faces)
        //if v_i was removed, add \frac{1}{deg(v_i)} \times residual(v_i) to the lifting value;
        // else add 0.
        //then divide this value by the numùber of these neighbours, and lift the position of the vertex.
        VectorInt lift = CGAL::NULL_VECTOR;
        int nNeighbours = vh->vertex_degree();

        hIt = h;
        do
        {
            Facet_handle f = hIt->facet();
            if(f->isSplittable())
                lift = lift + f->getResidual();
            hIt = hIt->next()->opposite();
        }
        while (hIt != h);

        PointInt quantPos = getQuantizedPos(vh->point());

        if (b_unlift)
            lift = -lift;

        PointInt newPosInt = quantPos + lift / nNeighbours / INV_GAMMA;

        for (unsigned i = 0; i < 3; ++i)
            assert(newPosInt[i] >= 0
                   && newPosInt[i] < 1 << (i_quantBits - i_curQuantizationId + LIFTING_NB_ADDITIONAL_BITS_GEOMETRY));

        vh->point() = getPos(newPosInt);
    }
}
