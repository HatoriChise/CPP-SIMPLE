#pragma once

#include "../config/simulationConfig.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

// 边界条件场：为后续离散方程组装提供按边、按索引的条件查询/写入
class BoundaryField {
public:
	BoundaryField(int ncx, int ncy, const BoudaryCondition* configs, std::size_t count);

	// 读写访问：竖边用 j 索引 (0..ncy-1)，横边用 i 索引 (0..ncx-1)
	BoudaryCondition& west(int j);   // X_MIN
	BoudaryCondition& east(int j);   // X_MAX
	BoudaryCondition& south(int i);  // Y_MIN
	BoudaryCondition& north(int i);  // Y_MAX

	const BoudaryCondition& west(int j) const;
	const BoudaryCondition& east(int j) const;
	const BoudaryCondition& south(int i) const;
	const BoudaryCondition& north(int i) const;

	int ncx() const { return ncx_; }
	int ncy() const { return ncy_; }

private:
	int ncx_{};
	int ncy_{};

	// 每条边分配一维数组，便于按索引快速取值
	std::vector<BoudaryCondition> west_;  // size = ncy
	std::vector<BoudaryCondition> east_;  // size = ncy
	std::vector<BoudaryCondition> south_; // size = ncx
	std::vector<BoudaryCondition> north_; // size = ncx

	// 从配置表提取单边默认条件
	BoudaryCondition default_west_{};
	BoudaryCondition default_east_{};
	BoudaryCondition default_south_{};
	BoudaryCondition default_north_{};

	void initialize_defaults(const BoudaryCondition* configs, std::size_t count);
	void allocate_edges();
};

