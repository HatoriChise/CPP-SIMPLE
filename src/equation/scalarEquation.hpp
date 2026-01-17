// src/euation/scalarEquation.hpp

#pragma once
#include "src/field/boundaryField.hpp"
#include "src/field/fluidProperties.hpp"
#include "src/field/scalarField.hpp"
#include "src/field/vectorField.hpp"
#include "src/grid/structuredMesh.hpp"

#include <boost/multi_array.hpp>

/**
 * @brief coefficient of the scalar equation
 * stores the coefficients for the discretized scalar transport equation
 */
struct COEF
{
    float aE;   // east coefficient
    float aW;   // west coefficient
    float aN;   // north coefficient
    float aS;   // south coefficient
    float aP;   // central coefficient
    float bsrc; // source term
};

/**
 * @brief scalar equation class
  Represents a scalar transport equation in the CFD solver.
  It holds references to the mesh, scalar field, velocity field,
  boundary conditions, and fluid properties needed to assemble and solve the equation.
  It also contains a coefficient matrix for the discretized equation.
 *
 */
class ScalarEquation
{
    // data for scalar equation
private:
    StructuredMesh &mesh_;
    ScalarField &scalarField_;
    VectorField &vectorField_;
    BoundaryField &boundaryField_;
    FluidPropertyField &fluidPropertyField_;

    boost::multi_array<COEF, 2> coefMatrix_; // coefficients array

public:
    // constructor
    ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField, VectorField &vectorField,
                   BoundaryField &boundaryField, FluidPropertyField &fluidPropertyField);

    // destructor
    ~ScalarEquation();

    // member functions

    /**
     * @brief initialize the coefficient matrix or clean it before assembly
     * 
     */
    void resetCoefficients();

    void addConvectionTerm();

    void addDiffusionTerm();

    void addSourceTerm();

    void applyBoundaries();

    void setRelaxation(float relaxationFactor);

    const boost::multi_array<COEF, 2>& getCoefMatrix() const
    {
        return coefMatrix_;
    }    
};