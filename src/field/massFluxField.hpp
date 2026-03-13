// src/field/massFluxField.hpp
#pragma once

#include <boost/multi_array.hpp>
#include "src/grid/structuredMesh.hpp"
#include "src/field/scalarField.hpp"
#include "src/field/vectorField.hpp"
#include "src/field/fluidProperties.hpp"

// Forward declaration
class ScalarEquation;

struct FaceFlux {
    float mE{0.0f}; // East mass flux
    float mW{0.0f}; // West mass flux
    float mN{0.0f}; // North mass flux
    float mS{0.0f}; // South mass flux
};

class MassFluxField {
public:
    MassFluxField(int ncx, int ncy);
    ~MassFluxField() = default;

    const FaceFlux& operator()(int i, int j) const { return data_[j][i]; }
    FaceFlux& operator()(int i, int j) { return data_[j][i]; }

    int ncx() const { return ncx_; }
    int ncy() const { return ncy_; }

    void updateFluxes(const StructuredMesh& mesh, 
                      const VectorField& vectorField, 
                      const FluidPropertyField& fluidPropertyField, 
                      const ScalarField* pressure,
                      const ScalarEquation* uEq, 
                      const ScalarEquation* vEq);

private:
    int ncx_;
    int ncy_;
    boost::multi_array<FaceFlux, 2> data_;

    enum class Face { East, West, North, South };

    float computeFaceMassFlux(int i, int j, Face face, 
                              const StructuredMesh& mesh,
                              const VectorField& vectorField,
                              const FluidPropertyField& fluidPropertyField,
                              const ScalarField* pressure,
                              const ScalarEquation* uEq,
                              const ScalarEquation* vEq) const;
};
