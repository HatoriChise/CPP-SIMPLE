// src/grid/structuredMesh.hpp
#pragma once

#include <Eigen/Dense>
#include <vector>
#include <stdexcept>
#include <fmt/core.h>
#include <boost/multi_array.hpp>   

#include "simulationConfig.hpp"

class StructuredMesh
{

    private:
    // cell size in each direction
    float dx_;
    float dy_;

    // cell centers in each direction
    std::vector<float> cellCentersX_;
    std::vector<float> cellCentersY_;

    // cell surfaces in each direction
    std::vector<float> nodeCoordinatesX_;
    std::vector<float> nodeCoordinatesY_;

    // for each cell, store the indices of its four faces
    boost::multi_array<struct CELL_FACE, 2> cellFaceID_; // shape: (numCellsY, numCellsX)



    public:
    // constructor
    StructuredMesh();
    // destructor
    ~StructuredMesh();

    // member functions
    void createVolumeMesh();
    void createBoundaryMesh();

    /**
     * @brief Get the Mesh Size object 对于均匀网格，网格尺寸相同
     * 
     * @return std::array<float, 3> 
     */
    std::array<float, 3> getMeshSize();

    void printMeshInfo(); // 打印网格信息


    inline const std::vector<float>& getCellCentersX() const { return cellCentersX_; }
    inline const std::vector<float>& getCellCentersY() const { return cellCentersY_; }
    
    inline const std::vector<float>& getNodeCoordinatesX() const { return nodeCoordinatesX_; }
    inline const std::vector<float>& getNodeCoordinatesY() const { return nodeCoordinatesY_; }

    inline const boost::multi_array<struct CELL_FACE, 2>& getCellFaceID() const { return cellFaceID_; }
};