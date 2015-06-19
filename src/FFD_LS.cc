#include "FFD_LS.h"

FFD_LS::FFD_LS(const Problem & problem, const std::string & directory, const Real & step, const Index & maxIterationsNo, const Real & tolerance, const bool & volume_constraint, const std::pair<Point, Point> & boundingBox, const std::pair<Index, Index> & sub, const Real & beta, const Real & armijoSlope)
    : FFD(problem, directory, step, maxIterationsNo, tolerance, volume_constraint, boundingBox, sub, armijoSlope), beta_(beta)
{
    //Scanning the border nodes
    Mesh::const_element_iterator       ref_el     = reference_mesh_.active_local_elements_begin();
    const Mesh::const_element_iterator ref_end_el = reference_mesh_.active_local_elements_end();
    
    for ( ; ref_el != ref_end_el; ++ref_el)
    {
        const Elem * elem = *ref_el;
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (elem->neighbor(side) == NULL)
            {
                border_ref_.push_back(*(elem->get_node(side)));
            }
        }
    }
    
    //Creating the B_ref Matrix
    Index K = CP_grid_.cols();
    Index L = CP_grid_.rows();
    Index LL = K * L;
    Index NB = border_ref_.size();
    
    B_x.resize(NB, LL);
    B_y.resize(NB, LL);
    
    for (Index i = 0; i < NB; ++i)
    {
        for (Index ll = 0; ll < LL; ++ll)
        {
            Index k = ll % K;
            Index l = ll / K;
            B_x(i, ll) = basisFunction(psi(border_ref_[i]), k, l) * (boundingBox_.second(0) - boundingBox_.first(0));
            B_y(i, ll) = basisFunction(psi(border_ref_[i]), k, l) * (boundingBox_.second(1) - boundingBox_.first(1));
        }
    }
    
    solver_x.compute( beta * B_x.transpose() * B_x + (1 - beta) * MatrixXd::Identity(LL, LL) );
    solver_y.compute( beta * B_y.transpose() * B_y + (1 - beta) * MatrixXd::Identity(LL, LL) );
}

void FFD_LS::computePerturbation(EquationSystems & perturbation, EquationSystems & stateAdj)
{
    Index K = CP_grid_.cols();
    Index L = CP_grid_.rows();
    Index LL = K * L;
    Index NB = border_ref_.size();
    
    VectorXd f1(NB);
    VectorXd f2(NB);
    
    Mesh::const_element_iterator       el     = mesh_->active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh_->active_local_elements_end();
    
    /*
     * Compute the average between the two adjacent edge normals for each boundary vertex.
     */
    std::map<Index, Point> normals;
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (elem->neighbor(side) == NULL)
            {
                Node * node           = elem->get_node(side);
                Node * following_node = elem->get_node((side + 1) % 3);
                
                Real deltax = (*following_node)(0) - (*node)(0);
                Real deltay = (*following_node)(1) - (*node)(1);
                
                // Somma vettoriale tra le due normali.
                normals[node->id()](0) += deltax;
                normals[node->id()](1) += deltay;
                
                normals[following_node->id()](0) += deltax;
                normals[following_node->id()](1) += deltay;
            }
        }
    }
    
    //
    
    el = mesh_->active_local_elements_begin();
    Index count = 0;
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (elem->neighbor(side) == NULL)
            {
                Node * node = elem->get_node(side);
                
                Real g = problem_.computeGradient(stateAdj, *node) + actual_lagrange_;
                
                Real nx = normals[node->id()](0);
                Real ny = normals[node->id()](1);
                
                Real modulus = std::sqrt(nx * nx + ny * ny);
                
                nx /= modulus;
                ny /= modulus;
                
                f1[count] = (*node)(0) - step_ * g * ny - border_ref_[count](0);
                f2[count] = (*node)(1) + step_ * g * nx - border_ref_[count](1);
                
                // std::cout << (*node)(0) << " " << (*node)(1) << " " << f1[count] << " " << f2[count] << ";" << std::endl;
                
                ++count;
            }
        }
    }
    
    VectorXd mu_x = solver_x.solve(B_x.transpose() * f1);
    VectorXd mu_y = solver_y.solve(B_y.transpose() * f2);
    
    for (Index ll = 0; ll < LL; ++ll)
    {
        Index k = ll % K;
        Index l = ll / K;
        gradJ_(L - 1 - l, k)(0) = (mu_(L - 1 - l, k)(0) - mu_x(ll)) / step_; // So that mu_ = - [mu_x, mu_y];
        gradJ_(L - 1 - l, k)(1) = (mu_(L - 1 - l, k)(1) - mu_y(ll)) / step_;
    }
}
