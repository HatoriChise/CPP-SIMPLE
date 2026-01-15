// src/grid/structuredGrid.cpp

#include "structuredMesh.hpp"

StructuredMesh::StructuredMesh()
{ 
    // constructor implementation (if needed)
    // compute cell sizes based on simulation config
    this->dx_ = Lx / ncx;
    this->dy_ = Ly / ncy;

    
}