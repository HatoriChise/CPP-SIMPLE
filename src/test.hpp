#pragma once

#include <fmt/core.h>
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>

#include "utils/formatter4eigen.h"

#include "config/simulationConfig.hpp"

#include "grid/structuredMesh.hpp"

#include "field/scalarField.hpp"
#include "field/vectorField.hpp"
#include "field/boundaryField.hpp"
#include "field/fluidProperties.hpp"
#include "equation/scalarEquation.hpp"

// define test function
void test();
void test_simple_algorithm();