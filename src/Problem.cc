#include "Problem.h"

Problem::Problem(Mesh mesh)
    : mesh_(std::make_shared<Mesh>(mesh)) {}
