#include "location_boundary.h"

bool isDirLike(PDEBoundaryType t) { return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented; }