/* C++ */

/**
 * @file   typedefs.h
 * @author Pasquale Claudio Africa <pasquale.africa@mail.polimi.it>, Luca Ratti <luca3.ratti@mail.polimi.it>, Abele Simona <abele.simona@mail.polimi.it>
 * @date   2015
 *
 * Questo file fa parte del progetto "ShapeOpt".
 *
 * @copyright Copyright © 2014 Pasquale Claudio Africa, Luca Ratti, Abele Simona. All rights reserved.
 * @copyright This project is released under the GNU General Public License.
 *
 * @brief Confronto tra alcune tecniche per l'ottimizzazione di forma.
 *
 */

#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <iostream>
#include <fstream>
#include <memory>

#include <boost/math/special_functions/binomial.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/distributed_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/vtk_io.h"
#include "libmesh/zero_function.h"

#include "GetPot.h"

using Real = double;    /**< @brief Typedef for real numbers. */
using Index = ptrdiff_t;    /**< @brief Typedef for indexing variables. */

using namespace boost::math;
using namespace Eigen;
using namespace libMesh;

using MatrixXp = Matrix<Point, Dynamic, Dynamic>;    /**< @brief Typedef for dense dynamic-sized matrices of points. */
using VectorXp = Matrix<Point, Dynamic, 1>      ;    /**< @brief Typedef for dense dynamic-sized column vectors of points. */
using MatrixXr = Matrix<Real, Dynamic, Dynamic> ;    /**< @brief Typedef for dense real-valued dynamic-sized matrices. */
using VectorXr = Matrix<Real, Dynamic,       1> ;    /**< @brief Typedef for dense real-valued dynamic-sized column vectors. */

using SparseXr = Eigen::SparseMatrix<Real>;    /**< @brief Typedef for sparse real-valued dynamic-sized matrices. */

#endif /* TYPEDEFS_H */
