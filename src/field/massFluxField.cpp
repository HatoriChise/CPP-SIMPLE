// src/field/massFluxField.cpp

#include <cmath>
#include <stdexcept>
#include "src/equation/scalarEquation.hpp"
#include "src/field/massFluxField.hpp"

MassFluxField::MassFluxField(int ncx, int ncy) : ncx_(ncx), ncy_(ncy)
{
    data_.resize(boost::extents[ncy_][ncx_]);
}

void MassFluxField::initalizeFluxes()
{
    for(int j = 0; j < ncy_; ++j)
    {
        for(int i = 0; i < ncx_; ++i)
        {
            data_[j][i].mE = 0.0f;
            data_[j][i].mW = 0.0f;
            data_[j][i].mN = 0.0f;
            data_[j][i].mS = 0.0f;
        }
    }
}
void MassFluxField::updateFluxes(const StructuredMesh &mesh, const VectorField &vectorField,
                                 const FluidPropertyField &fluidPropertyField,
                                 const BoundaryField &boundaryField, const ScalarField *pressure,
                                 const ScalarEquation *uEq, const ScalarEquation *vEq)
{
    for(int j = 0; j < ncy_; ++j)
    {
        for(int i = 0; i < ncx_; ++i)
        {
            data_[j][i].mE =
                computeFaceMassFlux(i, j, Face::East, mesh, vectorField, fluidPropertyField,
                                    boundaryField, pressure, uEq, vEq);
            data_[j][i].mW =
                computeFaceMassFlux(i, j, Face::West, mesh, vectorField, fluidPropertyField,
                                    boundaryField, pressure, uEq, vEq);
            data_[j][i].mN =
                computeFaceMassFlux(i, j, Face::North, mesh, vectorField, fluidPropertyField,
                                    boundaryField, pressure, uEq, vEq);
            data_[j][i].mS =
                computeFaceMassFlux(i, j, Face::South, mesh, vectorField, fluidPropertyField,
                                    boundaryField, pressure, uEq, vEq);
        }
    }
}

float MassFluxField::computeFaceMassFlux(int i, int j, Face face, const StructuredMesh &mesh,
                                         const VectorField &vectorField,
                                         const FluidPropertyField &fluidPropertyField,
                                         const BoundaryField &boundaryField,
                                         const ScalarField *pressure, const ScalarEquation *uEq,
                                         const ScalarEquation *vEq) const
{
    auto meshSize = mesh.getMeshSize();
    const auto& bc = boundaryField;
    float dx = meshSize[0];
    float dy = meshSize[1];
    float volume = dx * dy;
    float rho = fluidPropertyField(i, j).rho;

    float linearVelocity = 0.0f;

    // 先判断是否在边界，如果在边界，根据边界条件调整线速度计算方式。然后return
    if (i == 0 && face == Face::West) // XMIN处
    {   
        if (bc.west(j).velocityType == INLET)
        {
            return rho * bc.west(j).VelocityValue[0] * dy;
        }
        else if (bc.west(j).velocityType == OUTLET)
        {
            return rho * vectorField.u()(i, j) * dy; // 使用当前单元的速度
        }
        else if (bc.west(j).velocityType == WALL)
        {
            return 0.0f; // 壁面无流动
        }
    }

    if (i == ncx_ - 1 && face == Face::East) // XMAX处
    {
        if (bc.east(j).velocityType == INLET)
        {
            return rho * bc.east(j).VelocityValue[0] * dy;
        }
        else if (bc.east(j).velocityType == OUTLET)
        {
            return rho * vectorField.u()(i, j) * dy; // 使用当前单元的速度
        }
        else if (bc.east(j).velocityType == WALL)
        {
            return 0.0f; // 壁面无流动
        }
    }

    if (j == 0 && face == Face::South) // YMIN处
    {
        if (bc.south(i).velocityType == INLET)
        {
            return rho * bc.south(i).VelocityValue[1] * dx;
        }
        else if (bc.south(i).velocityType == OUTLET)
        {
            return rho * vectorField.v()(i, j) * dx; // 使用当前单元的速度
        }
        else if (bc.south(i).velocityType == WALL)
        {
            return 0.0f; // 壁面无流动
        }
    }

    if (j == ncy_ - 1 && face == Face::North) // YMAX处
    {
        if (bc.north(i).velocityType == INLET)
        {
            return rho * bc.north(i).VelocityValue[1] * dx;
        }
        else if (bc.north(i).velocityType == OUTLET)
        {
            return rho * vectorField.v()(i, j) * dx; // 使用当前单元
        }
        else if (bc.north(i).velocityType == WALL)
        {
            return 0.0f; // 壁面无流动
        }
    }


    float area = (face == Face::East || face == Face::West) ? dy : dx;
    switch(face)
    {
    case Face::East:
        linearVelocity = 0.5f * (vectorField.u()(i, j) + vectorField.u()(i + 1, j));
        break;
    case Face::West:
        linearVelocity = 0.5f * (vectorField.u()(i - 1, j) + vectorField.u()(i, j));
        break;
    case Face::North:
        linearVelocity = 0.5f * (vectorField.v()(i, j) + vectorField.v()(i, j + 1));
        break;
    case Face::South:
        linearVelocity = 0.5f * (vectorField.v()(i, j - 1) + vectorField.v()(i, j));
        break;
    }

    if(pressure == nullptr || uEq == nullptr || vEq == nullptr)
    {
        return rho * linearVelocity * area;
    }

    const boost::multi_array<COEF, 2> &coefU = uEq->getCoefMatrix();
    const boost::multi_array<COEF, 2> &coefV = vEq->getCoefMatrix();

    float aP = 0.0f;
    float aP_nabor = 0.0f;

    switch(face)
    {
    case Face::East:
        aP = coefU[j][i].aP;
        aP_nabor = (i + 1 < ncx_) ? coefU[j][i + 1].aP : aP; // 边界用aP
        break;
    case Face::West:
        aP = coefU[j][i].aP;
        aP_nabor = (i > 0) ? coefU[j][i - 1].aP : aP; // 边界用aP
        break;
    case Face::North:
        aP = coefV[j][i].aP;
        aP_nabor = (j + 1 < ncy_) ? coefV[j + 1][i].aP : aP; // 边界用aP
        break;
    case Face::South:
        aP = coefV[j][i].aP;
        aP_nabor = (j > 0) ? coefV[j - 1][i].aP : aP; // 边界用aP ✅
        break;
    }

    if(std::fabs(aP) < 1e-10f || std::fabs(aP_nabor) < 1e-10f)
    {
        return rho * linearVelocity * area;
    }

    float d_P = volume / aP;
    float d_nabor = volume / aP_nabor;
    float d_face = 0.5f * (d_P + d_nabor);

    // gradP_P = (P_W - P_E) / 2*dx 或 (P_S - P_N) / 2*dy
    float gradP_P = 0.0f;
    switch(face)
    {
    case Face::East:
    case Face::West:
        if(i > 0 && i < ncx_ - 1)
        {
            gradP_P =
                (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);
        }
        else if(i == 0)
        {
            gradP_P = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else
        {
            gradP_P = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        break;
    case Face::North:
    case Face::South:
        if(j > 0 && j < ncy_ - 1)
        {
            gradP_P =
                (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);
        }
        else if(j == 0)
        {
            gradP_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else
        {
            gradP_P = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        break;
    }

    // gradP_nabor = (P_EE - P_P) / 2*dx 或 (P_NN - P_P) / 2*dy
    // (P_WW - P_P) / 2*dx 或 (P_SS - P_P) / 2*dy
    float gradP_nabor = 0.0f;
    switch(face)
    {
    case Face::East:
        if(i < ncx_ - 2)
        {
            gradP_nabor =
                (pressure->operator()(i + 2, j) - pressure->operator()(i, j)) / (2.0f * dx);
        }
        else if(i == ncx_ - 2)
        {
            gradP_nabor = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        break;
    case Face::West:
        if(i > 1)
        {
            gradP_nabor =
                (pressure->operator()(i, j) - pressure->operator()(i - 2, j)) / (2.0f * dx);
        }
        else
        {
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        break;
    case Face::North:
        if(j < ncy_ - 2)
        {
            gradP_nabor =
                (pressure->operator()(i, j + 2) - pressure->operator()(i, j)) / (2.0f * dy);
        }
        else
        {
            gradP_nabor = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        break;
    case Face::South:
        if(j > 1)
        {
            gradP_nabor =
                (pressure->operator()(i, j) - pressure->operator()(i, j - 2)) / (2.0f * dy);
        }
        else
        {
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        break;
    }

    // 面上的压力梯度
    float gradP_face = 0.0f;
    switch(face)
    {
    case Face::East:
        if(i + 1 < ncx_)
        {
            gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else if(i - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::West:
        if(i - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        else if(i + 1 < ncx_)
        {
            gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::North:
        if(j + 1 < ncy_)
        {
            gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else if(j - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::South:
        if(j - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        else if(j + 1 < ncy_)
        {
            gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    }

    float correctedVelocity =
        linearVelocity + 0.5f * d_P * gradP_P + 0.5f * d_nabor * gradP_nabor - d_face * gradP_face;

    // Debug output (formatted)
    std::cout << "===========================================================" << std::endl;
    std::cout << "Debug Info for cell (" << i << ", " << j << ") - Face: ";
    switch(face)
    {
    case Face::East:
        std::cout << "East";
        break;
    case Face::West:
        std::cout << "West";
        break;
    case Face::North:
        std::cout << "North";
        break;
    case Face::South:
        std::cout << "South";
        break;
    }
    std::cout << "\nDEBUG_nabor:" << " volume:" << volume << " aP_nabor: " << aP_nabor
              << " d_nabor: " << d_nabor << std::endl;
    std::cout << "\n"
              << " d_P       : " << d_P << " gradP_P   : " << gradP_P
              << " d_p * gradP_P: " << d_P * gradP_P << "\n"
              << " d_nabor   : " << d_nabor << " gradP_nabor: " << gradP_nabor
              << " d_nabor * gradP_nabor: " << d_nabor * gradP_nabor << "\n"
              << " d_face    : " << d_face << " gradP_face : " << gradP_face
              << " d_face * gradP_face: " << d_face * gradP_face << "\n"
              << " linearVelocity: " << linearVelocity << "\n"
              << " correctedVelocity: " << correctedVelocity << "\n"
              << std::endl;

    return rho * correctedVelocity * area;
}
