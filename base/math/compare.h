#pragma once

// https://www.learncpp.com/cpp-tutorial/relational-operators-and-floating-point-comparisons/
// Return true if the difference between a and b is within epsilon percent of the larger of a and b
bool approximatelyEqualRel(double a, double b, double relEpsilon);
// Return true if the difference between a and b is less than or equal to absEpsilon, or within relEpsilon percent of
// the larger of a and b
bool approximatelyEqualAbsRel(double a, double b, double absEpsilon, double relEpsilon);