#pragma once

#include "base/pch.h"

#include <random>

void add_random_number(double* field, int length, double min, double max, int seed = -1);
void add_random_number(field2& field, double min, double max, int seed = -1);
void add_random_number(field3& field, double min, double max, int seed = -1);