// src/euation/scalarEquation.hpp

#pragma once
#include "src/field/boundaryField.hpp"
#include "src/field/fluidProperties.hpp"
#include "src/field/scalarField.hpp"
#include "src/field/vectorField.hpp"
#include "src/grid/structuredMesh.hpp"

#include <boost/multi_array.hpp>
#include <iostream>

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

    int direction_;  // 0=u动量, 1=v动量, -1=标量方程(默认)

    // 计算界面质量通量（当前使用线性插值，后续可升级为Rhie-Chow）
    // face: 0=东, 1=西, 2=北, 3=南
    float computeFaceMassFlux(int i, int j, int face, const ScalarField* pressure = nullptr) const;

    // 计算界面扩散系数（调和平均）
    // mu_owner: owner 单元的粘度
    // mu_neighbor: neighbor 单元的粘度
    // distance: owner 到 neighbor 的中心距离
    // area: 界面面积
    float computeFaceDiffusionCoefficient(float mu_owner, float mu_neighbor, 
                                           float distance, float area) const;

public:
    // constructor
    // direction: 0=u动量方程, 1=v动量方程, -1=标量方程(默认)
    ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField, VectorField &vectorField,
                   BoundaryField &boundaryField, FluidPropertyField &fluidPropertyField,
                   int direction = -1);

    // destructor
    ~ScalarEquation();

    // member functions

    /**
     * @brief initialize the coefficient matrix or clean it before assembly
     * 
     */
    void resetCoefficients();

    void addConvectionTerm(const ScalarField* pressure = nullptr);

    void addDiffusionTerm();

    void addSourceTerm();

    void applyBoundaries();

    // 添加压力梯度源项，需要传入压力场
    void addPressureGradient(const ScalarField& pressure);

    void setRelaxation(float relaxationFactor);

    const boost::multi_array<COEF, 2>& getCoefMatrix() const
    {
        return coefMatrix_;
    }

    // 保存系数到文件
    void saveCoefficientsToFile(const std::string& filename) const;
};