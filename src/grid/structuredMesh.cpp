// src/grid/structuredMesh.cpp

#include "structuredMesh.hpp"

#include <filesystem>
#include <fstream>

StructuredMesh::StructuredMesh()
{
    if (ncx <= 0 || ncy <= 0)
    {
        throw std::invalid_argument("ncx and ncy must be positive");
    }

    dx_ = Lx / static_cast<float>(ncx);
    dy_ = Ly / static_cast<float>(ncy);

    createVolumeMesh();
    createBoundaryMesh();
}

StructuredMesh::~StructuredMesh() = default;

void StructuredMesh::createVolumeMesh()
{
    cellCentersX_.resize(ncx);
    for (int i = 0; i < ncx; ++i)
    {
        cellCentersX_[static_cast<std::size_t>(i)] = (static_cast<float>(i) + 0.5f) * dx_;
    }

    cellCentersY_.resize(ncy);
    for (int j = 0; j < ncy; ++j)
    {
        cellCentersY_[static_cast<std::size_t>(j)] = (static_cast<float>(j) + 0.5f) * dy_;
    }

    nodeCoordinatesX_.resize(ncx + 1);
    for (int i = 0; i <= ncx; ++i)
    {
        nodeCoordinatesX_[static_cast<std::size_t>(i)] = static_cast<float>(i) * dx_;
    }

    nodeCoordinatesY_.resize(ncy + 1);
    for (int j = 0; j <= ncy; ++j)
    {
        nodeCoordinatesY_[static_cast<std::size_t>(j)] = static_cast<float>(j) * dy_;
    }

    cellFaceID_.resize(boost::extents[ncy][ncx]);
}

void StructuredMesh::createBoundaryMesh()
{
    if (cellFaceID_.num_elements() == 0)
    {
        cellFaceID_.resize(boost::extents[ncy][ncx]);
    }

    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            CELL_FACE faces{INTERIOR, INTERIOR, INTERIOR, INTERIOR};

            if (i == 0)
            {
                faces.WEST = X_MIN;
            }
            if (i == ncx - 1)
            {
                faces.EAST = X_MAX;
            }
            if (j == 0)
            {
                faces.SOUTH = Y_MIN;
            }
            if (j == ncy - 1)
            {
                faces.NORTH = Y_MAX;
            }

            cellFaceID_[static_cast<std::size_t>(j)][static_cast<std::size_t>(i)] = faces;
        }
    }
}

std::array<float, 3> StructuredMesh::getMeshSize()
{
    return {dx_, dy_, dx_ * dy_};
}

void StructuredMesh::saveMeshInfo()
{
    namespace fs = std::filesystem;

    fs::path outputDir{"data"};
    if (!fs::exists(outputDir))
    {
        outputDir = fs::path("..") / "data";
    }
    fs::create_directories(outputDir);

    fs::path outputFile = outputDir / "mesh_info.json";
    std::ofstream out(outputFile);
    if (!out)
    {
        throw std::runtime_error(fmt::format("Failed to open {} for writing", outputFile.string()));
    }

    out << "{\n";
    out << fmt::format("  \"ncx\": {},\n", ncx);
    out << fmt::format("  \"ncy\": {},\n", ncy);
    out << fmt::format("  \"Lx\": {:.8f},\n", Lx);
    out << fmt::format("  \"Ly\": {:.8f},\n", Ly);
    out << fmt::format("  \"dx\": {:.8f},\n", dx_);
    out << fmt::format("  \"dy\": {:.8f},\n", dy_);
    out << fmt::format("  \"cell_area\": {:.8f},\n", dx_ * dy_);

    auto dumpVector = [&out](const std::vector<float>& values, const std::string& name, bool add_comma) {
        out << fmt::format("  \"{}\": [", name);
        for (std::size_t idx = 0; idx < values.size(); ++idx)
        {
            if (idx != 0)
            {
                out << ", ";
            }
            out << fmt::format("{:.8f}", values[idx]);
        }
        out << "]" << (add_comma ? ",\n" : "\n");
    };

    dumpVector(cellCentersX_, "cell_centers_x", true);
    dumpVector(cellCentersY_, "cell_centers_y", true);
    dumpVector(nodeCoordinatesX_, "node_x", true);
    dumpVector(nodeCoordinatesY_, "node_y", false);

    out << "}\n";
}