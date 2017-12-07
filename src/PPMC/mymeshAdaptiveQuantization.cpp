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


/**
  * Determine the value of the Kg coeficient according to the method of Ho Lee
  * in the paper:
  * "Rate-distortion optimization for progressive compression of 3D mesh with color attributes"
  */
float MyMesh::determineKg()
{
    float f_surface = 0;
    for (Facet_iterator fit = facets_begin(); fit != facets_end(); ++fit)
        f_surface += faceSurface(fit->facet_begin());

    float f_Kg = f_bbVolume * f_adaptQuantRescaling / (f_surface * size_of_vertices());

    return f_Kg;
}


/**
  * Determine the symbol of each subdivision cell.
  */
std::map<unsigned, unsigned> MyMesh::determineCellSymbols(Halfedge_handle heh_v, bool b_compression)
{
#if 1
    float d[3][2] = {0};
    Point pos = heh_v->vertex()->point();

    // Compute the distances.
    Halfedge_const_handle heh_init = heh_v->opposite(), heh = heh_init;
    do
    {
        Vertex_const_handle vh_neighbor = heh->vertex();
        Point neightborPos =  b_compression && vh_neighbor->isConquered() ?
                              vh_neighbor->getOldPos() : vh_neighbor->point();

        unsigned i_level = i_quantBits - i_curQuantizationId;
        if (vh_neighbor->isConquered())
            i_level++;

        for (unsigned i = 0; i < 3; ++i)
        {
            float f_dist = sqrt((neightborPos - pos).squared_length());
            if (neightborPos[i] <= pos[i])
                d[i][0] += i_level * f_dist;
            else
                d[i][1] += i_level * f_dist;
        }
        heh = heh->prev()->opposite();
    }
    while (heh != heh_init);

    // Compute the weigths of unbalancing.
    std::multimap<float, unsigned> axesOrder;
    for (unsigned i = 0; i < 3; ++i)
        axesOrder.insert(std::pair<float, unsigned>(fabs(d[i][0] / (d[i][0] + d[i][1]) - 0.5), i));

    unsigned axesWeigth[3];
    unsigned i_weight = 1;
    for (std::multimap<float, unsigned>::iterator it = axesOrder.begin(); it != axesOrder.end(); it++)
        axesWeigth[it->second] = i_weight++;

    // Compute the priority of each cell.
    std::multimap<float, unsigned> cellOrder;
    for (unsigned i = 0; i < 8; ++i)
    {
        float f_priority = 0;
        unsigned cellcoresp[][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0},
                                    {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}};
        for (unsigned j = 0; j < 3; ++j)
            f_priority += axesWeigth[j] * d[j][cellcoresp[i][j]];

        cellOrder.insert(std::pair<float, unsigned>(f_priority, i));
    }

    // Determine the mapping between the cell number and the new symbol.
    std::map<unsigned, unsigned> cellMap;
    {
        unsigned i = 0;
        for (std::multimap<float, unsigned>::iterator it = cellOrder.begin();
             it != cellOrder.end(); it++, i++)
        {
            if (b_compression)
                cellMap.insert(std::pair<unsigned, unsigned>(it->second, 7 - i));
            else
                cellMap.insert(std::pair<unsigned, unsigned>(7 - i, it->second));
        }
    }

#else

    Point pos = heh_v->vertex()->point();

    std::vector<Vertex_const_handle> polygon;
    polygon.reserve(heh_v->vertex_degree());

    // Compute the distances.
    Halfedge_const_handle heh_init = heh_v->opposite(), heh = heh_init;
    do
        polygon.push_back(heh->vertex());
    while (heh != heh_init);

    Point b = barycenter(polygon);

    Vector cellPos[] = {Vector(-1, -1, -1), Vector(1, -1, -1), Vector(-1, 1, -1), Vector(1, 1, -1),
                        Vector(-1, -1, 1),  Vector(1, -1, 1),  Vector(-1, 1, 1),  Vector(1, 1, 1)};

    std::multimap<float, unsigned> cellOrder;
    for (unsigned i = 0; i < 8; ++i)
        cellOrder.insert(std::pair<float, unsigned>(
            sqrt((b - pos + cellPos[i] * f_quantStep * (1 << i_curQuantizationId - 1)).squared_length()),
            i));

    std::map<unsigned, unsigned> cellMap;
    {
        unsigned i = 0;
        for (std::multimap<float, unsigned>::iterator it = cellOrder.begin();
             it != cellOrder.end(); it++, i++)
        {
            if (b_compression)
                cellMap.insert(std::pair<unsigned, unsigned>(it->second, 7 - i));
            else
                cellMap.insert(std::pair<unsigned, unsigned>(7 - i, it->second));
        }
    }
#endif

    return cellMap;
}
