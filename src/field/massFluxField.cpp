// src/field/massFluxField.cpp

#include "src/field/massFluxField.hpp"
#include "src/equation/scalarEquation.hpp"
#include <cmath>
#include <stdexcept>

MassFluxField::MassFluxField(int ncx, int ncy) : ncx_(ncx), ncy_(ncy) {
    data_.resize(boost::extents[ncy_][ncx_]);
}

void MassFluxField::updateFluxes(const StructuredMesh& mesh, 
                                 const VectorField& vectorField, 
                                 const FluidPropertyField& fluidPropertyField, 
                                 const ScalarField* pressure,
                                 const ScalarEquation* uEq, 
                                 const ScalarEquation* vEq) {
    for(int j = 0; j < ncy_; ++j) {
        for(int i = 0; i < ncx_; ++i) {
            data_[j][i].mE = computeFaceMassFlux(i, j, Face::East, mesh, vectorField, fluidPropertyField, pressure, uEq, vEq);
            data_[j][i].mW = computeFaceMassFlux(i, j, Face::West, mesh, vectorField, fluidPropertyField, pressure, uEq, vEq);
            data_[j][i].mN = computeFaceMassFlux(i, j, Face::North, mesh, vectorField, fluidPropertyField, pressure, uEq, vEq);
            data_[j][i].mS = computeFaceMassFlux(i, j, Face::South, mesh, vectorField, fluidPropertyField, pressure, uEq, vEq);
        }
    }
}

float MassFluxField::computeFaceMassFlux(int i, int j, Face face, 
                                         const StructuredMesh& mesh,
                                         const VectorField& vectorField,
                                         const FluidPropertyField& fluidPropertyField,
                                         const ScalarField* pressure,
                                         const ScalarEquation* uEq,
                                         const ScalarEquation* vEq) const {
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0]; 
    float dy = meshSize[1]; 
    float volume = dx * dy; 
    float rho = fluidPropertyField(i, j).rho;

    float linearVelocity = 0.0f;
    float area = (face == Face::East || face == Face::West) ? dy : dx;
    switch(face) {
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

    if(pressure == nullptr || uEq == nullptr || vEq == nullptr) {
        return rho * linearVelocity * area;
    }

    const boost::multi_array<COEF, 2>& coefU = uEq->getCoefMatrix();
    const boost::multi_array<COEF, 2>& coefV = vEq->getCoefMatrix();

    float aP = 0.0f;
    float aP_nabor = 0.0f;
    
    switch(face) {
        case Face::East: 
            aP = coefU[j][i].aP;
            aP_nabor = coefU[j][i + 1].aP;
            break;
        case Face::West: 
            aP = coefU[j][i].aP;
            aP_nabor = coefU[j][i - 1].aP;
            break;
        case Face::North: 
            aP = coefV[j][i].aP;
            aP_nabor = coefV[j + 1][i].aP;
            break;
        case Face::South: 
            aP = coefV[j][i].aP;
            aP_nabor = coefV[j - 1][i].aP;
            break;
    }

    if(std::fabs(aP) < 1e-10f || std::fabs(aP_nabor) < 1e-10f) {
        return rho * linearVelocity * area;
    }

    float d_P = volume / aP;
    float d_nabor = volume / aP_nabor;
    float d_face = 0.5f * (d_P + d_nabor);

    float gradP_P = 0.0f;
    switch(face) {
        case Face::East: 
        case Face::West: 
            if(i > 0 && i < ncx_ - 1) {
                gradP_P = (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);
            } else if(i == 0) {
                gradP_P = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
            } else {
                gradP_P = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            }
            break;
        case Face::North: 
        case Face::South: 
            if(j > 0 && j < ncy_ - 1) {
                gradP_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);
            } else if(j == 0) {
                gradP_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
            } else {
                gradP_P = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            }
            break;
    }
    
    float gradP_nabor = 0.0f;
    switch(face) {
        case Face::East: 
            if(i < ncx_ - 2) {
                gradP_nabor = (pressure->operator()(i + 2, j) - pressure->operator()(i, j)) / (2.0f * dx);
            } else {
                gradP_nabor = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
            }
            break;
        case Face::West: 
            if(i > 1) {
                gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i - 2, j)) / (2.0f * dx);
            } else {
                gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            }
            break;
        case Face::North: 
            if(j < ncy_ - 2) {
                gradP_nabor = (pressure->operator()(i, j + 2) - pressure->operator()(i, j)) / (2.0f * dy);
            } else {
                gradP_nabor = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
            }
            break;
        case Face::South: 
            if(j > 1) {
                gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i, j - 2)) / (2.0f * dy);
            } else {
                gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            }
            break;
    }
    
    float gradP_face = 0.0f;
    switch(face) {
        case Face::East:
            if(i + 1 < ncx_) {
                gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
            } else if(i - 1 >= 0) {
                gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            } else {
                gradP_face = 0.0f;
            }
            break;
        case Face::West:
            if(i - 1 >= 0) {
                gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            } else if(i + 1 < ncx_) {
                gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
            } else {
                gradP_face = 0.0f;
            }
            break;
        case Face::North:
            if(j + 1 < ncy_) {
                gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
            } else if(j - 1 >= 0) {
                gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            } else {
                gradP_face = 0.0f;
            }
            break;
        case Face::South:
            if(j - 1 >= 0) {
                gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            } else if(j + 1 < ncy_) {
                gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
            } else {
                gradP_face = 0.0f;
            }
            break;
    }

    float correctedVelocity = linearVelocity + 0.5f * d_P * gradP_P + 0.5f * d_nabor * gradP_nabor - d_face * gradP_face;

    return rho * correctedVelocity * area;
}
