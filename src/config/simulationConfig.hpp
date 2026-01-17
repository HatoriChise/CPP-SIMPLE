// src/config/SimulationConfig.h
#pragma once
#include <array>
#include <string>

enum FACE_POSITION
{
    INTERIOR,
    X_MIN,
    X_MAX,
    Y_MIN,
    Y_MAX
};
enum BOUNDARY_TYPE_U
{
    NONE,
    WALL,
    INLET,
    OUTLET
};
enum BOUNDARY_TYPE_T
{
    DIRICHLET,
    NEUMANN,
    ROBIN
};

struct CELL_FACE
{
    FACE_POSITION EAST, WEST, NORTH, SOUTH;
};

struct BoudaryCondition
{
    FACE_POSITION position;
    BOUNDARY_TYPE_U velocityType;
    std::array<float, 2> VelocityValue; // u, v
    BOUNDARY_TYPE_T temperatureType;
    float TemperatureValue;
    float heatFluxValue;
};

// === 几何与网格 ===
constexpr float Lx = 1.0; // 腔体宽度
constexpr float Ly = 1.0; // 腔体高度
constexpr int ncx = 40;   // cell number in all directions
constexpr int ncy = 40;


// === 边界条件 ===
constexpr float lid_velocity = 1.0; // 顶盖速度（U

// === 物理参数 ===
constexpr float Re = 100.0;                  // 雷诺数
constexpr float rho = 1.0;                   // 密度（通常设为1）
constexpr float nu = lid_velocity * Lx / Re; // 运动粘度 = U*L / Re

// === 求解控制 ===
constexpr int max_iterations = 10000;
constexpr float continuity_tolerance = 1e-6;
constexpr float momentum_tolerance = 1e-6;

// === 输出控制 ===
const std::string output_filename = "results.dat";
constexpr int output_interval = 100; // 每隔多少步输出一次结果


// === 边界条件 ===
const BoudaryCondition boundaryInfo[4] = 
{
    {X_MIN, WALL, {0.0, 0.0}, DIRICHLET, 273.0, 0.0},
    {X_MAX, WALL, {0.0, 0.0}, DIRICHLET, 273.0, 0.0},
    {Y_MIN, WALL, {0.0, 0.0}, DIRICHLET, 273.0, 0.0},
    {Y_MAX, WALL, {lid_velocity, 0.0}, DIRICHLET, 273.0, 0.0}};