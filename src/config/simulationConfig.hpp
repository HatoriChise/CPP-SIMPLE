// src/config/SimulationConfig.h
#pragma once
#include <string>

enum FACE_POSITION {INTERIOR, X_MIN, X_MAX, Y_MIN, Y_MAX};
enum BOUNDARY_TYPE_U {NONE, WALL, INLET, OUTLET};
enum BOUNDARY_TYPE_T {DIRICHLET, NEUMANN, ROBIN};

struct CELL_FACE
{
    FACE_POSITION EAST, WEST, NORTH, SOUTH;
};

// === 几何与网格 ===
float Lx = 1.0; // 腔体宽度
float Ly = 1.0; // 腔体高度
int ncx = 40;     // cell number in all directions
int ncy = 40;   

// === 物理参数 ===
float Re = 100.0; // 雷诺数
float rho = 1.0;  // 密度（通常设为1）
float nu;         // 运动粘度 = U*L / Re

// === 边界条件 ===
float lid_velocity = 1.0; // 顶盖速度（U）

// === 求解控制 ===
int max_iterations = 10000;
float continuity_tolerance = 1e-6;
float momentum_tolerance = 1e-6;

// === 输出控制 ===
std::string output_filename = "results.dat";
int output_interval = 100;