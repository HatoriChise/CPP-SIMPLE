#include "boundaryField.hpp"

namespace {

// 缺省边界，用于未提供配置时占位
constexpr BoudaryCondition kEmptyBoundary{
	INTERIOR, NONE, {0.0f, 0.0f}, DIRICHLET, 0.0f, 0.0f};

}  // namespace

BoundaryField::BoundaryField(int ncx, int ncy, const BoudaryCondition* configs, std::size_t count)
	: ncx_(ncx), ncy_(ncy) {
	assert(ncx_ > 0 && ncy_ > 0);
	initialize_defaults(configs, count);
	allocate_edges();
}

void BoundaryField::initialize_defaults(const BoudaryCondition* configs, std::size_t count) {
	// 先设为空，再用配置覆盖
	default_west_ = default_east_ = default_south_ = default_north_ = kEmptyBoundary;

	if (!configs || count == 0) return;

	for (std::size_t idx = 0; idx < count; ++idx) {
		const auto& bc = configs[idx];
		switch (bc.position) {
		case X_MIN: default_west_ = bc; break;
		case X_MAX: default_east_ = bc; break;
		case Y_MIN: default_south_ = bc; break;
		case Y_MAX: default_north_ = bc; break;
		default: break;
		}
	}
}

void BoundaryField::allocate_edges() {
	west_.assign(ncy_, default_west_);
	east_.assign(ncy_, default_east_);
	south_.assign(ncx_, default_south_);
	north_.assign(ncx_, default_north_);
}

BoudaryCondition& BoundaryField::west(int j) {
	assert(j >= 0 && j < ncy_);
	return west_[static_cast<std::size_t>(j)];
}

BoudaryCondition& BoundaryField::east(int j) {
	assert(j >= 0 && j < ncy_);
	return east_[static_cast<std::size_t>(j)];
}

BoudaryCondition& BoundaryField::south(int i) {
	assert(i >= 0 && i < ncx_);
	return south_[static_cast<std::size_t>(i)];
}

BoudaryCondition& BoundaryField::north(int i) {
	assert(i >= 0 && i < ncx_);
	return north_[static_cast<std::size_t>(i)];
}

const BoudaryCondition& BoundaryField::west(int j) const {
	assert(j >= 0 && j < ncy_);
	return west_[static_cast<std::size_t>(j)];
}

const BoudaryCondition& BoundaryField::east(int j) const {
	assert(j >= 0 && j < ncy_);
	return east_[static_cast<std::size_t>(j)];
}

const BoudaryCondition& BoundaryField::south(int i) const {
	assert(i >= 0 && i < ncx_);
	return south_[static_cast<std::size_t>(i)];
}

const BoudaryCondition& BoundaryField::north(int i) const {
	assert(i >= 0 && i < ncx_);
	return north_[static_cast<std::size_t>(i)];
}

