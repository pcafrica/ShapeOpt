#include "FFD.h"

FFD::FFD(const Problem & problem, const std::string & directory, const Real & step, const Index & maxIterationsNo, const Real & tolerance, const bool & volume_constraint, const std::pair<Point, Point> & boundingBox, const std::pair<Index, Index> & sub, const Real & armijoSlope)
    : ShapeOptimization(problem, directory, step, maxIterationsNo, tolerance, volume_constraint, armijoSlope), reference_mesh_(*problem_.get_mesh()), boundingBox_(boundingBox), sub_(sub), firstTime_(true)
{
    // Assemble the control points grid.
    CP_grid_.resize(sub_.second + 1, sub_.first + 1);
    
    Real xIncrement = (boundingBox_.second(0) - boundingBox_.first(0)) / sub_.first;
    Real yIncrement = (boundingBox_.second(1) - boundingBox_.first(1)) / sub_.second;
    
    for ( Index i = 0; i < CP_grid_.rows(); ++i )
    {
        for ( Index j = 0; j < CP_grid_.cols(); ++j )
        {
            CP_grid_(CP_grid_.rows() - i - 1, j)(0) = boundingBox_.first(0) + xIncrement * ((int) j);
            CP_grid_(CP_grid_.rows() - i - 1, j)(1) = boundingBox_.first(1) + yIncrement * ((int) i);
            CP_grid_(CP_grid_.rows() - i - 1, j)(2) = 0.0;
        }
    }
    
    // Update the bounding box.
    boundingBox_.first  = Point(CP_grid_(CP_grid_.rows() - 1, 0)(0), CP_grid_(CP_grid_.rows() - 1, 0)(1), CP_grid_(CP_grid_.rows() - 1, 0)(2));
    boundingBox_.second = Point(CP_grid_(0, CP_grid_.cols() - 1)(0), CP_grid_(0, CP_grid_.cols() - 1)(1), CP_grid_(0, CP_grid_.cols() - 1)(2));
    
    // Initialize mu_ and gradJ_.
    mu_    = MatrixXp::Zero(CP_grid_.rows(), CP_grid_.cols());
    gradJ_ = mu_;
}

void FFD::computePerturbation(EquationSystems & perturbation, EquationSystems & stateAdj)
{
    const unsigned int dim = reference_mesh_.mesh_dimension();
    LinearImplicitSystem & system = stateAdj.get_system<LinearImplicitSystem>(problem_.get_name());
    
    const DofMap& dof_map = system.get_dof_map();
    
    FEType fe_type = dof_map.variable_type(0);
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    Index quadNodesNo = qface.n_points();
    Index count = 0;
    
    if ( firstTime_ )
    {
        firstTime_ = false;
        
        Mesh::const_element_iterator       ref_el     = reference_mesh_.active_local_elements_begin();
        const Mesh::const_element_iterator ref_end_el = reference_mesh_.active_local_elements_end();
        
        // Conta numero di lati di bordo e numero di nodi di quadratura per ogni lato.
        Index edgesNo = 0;
        
        for ( ; ref_el != ref_end_el; ++ref_el)
        {
            const Elem * elem = *ref_el;
            
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                if (elem->neighbor(side) == NULL)
                {
                    ++edgesNo;
                }
            }
        }
        
        // Nodi di quadratura nella mesh di riferimento.
        reference_nodes_.resize( edgesNo * quadNodesNo );
        
        ref_el = reference_mesh_.active_local_elements_begin(); // Resetta iteratore.
        
        for ( ; ref_el != ref_end_el; ++ref_el)
        {
            const Elem * elem = *ref_el;
            
            for (unsigned int side = 0; side < elem->n_sides(); side++)
            {
                if (elem->neighbor(side) == NULL)
                {
                    const std::vector<Point> & qface_point = fe_face->get_xyz();
                    fe_face->reinit(elem, side);
                    
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        reference_nodes_(count * quadNodesNo + qp) = qface_point[qp];
                    }
                    
                    ++count;
                }
            }
        }
    }
    
    // Gradiente del funzionale costo.
    Mesh::const_element_iterator       el     = mesh_->active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh_->active_local_elements_end();
    
    count = 0;
    
    for ( ; el != end_el; ++el)
    {
        const Elem * elem = *el;
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
        {
            if (elem->neighbor(side) == NULL)
            {
                const std::vector<Real> & JxW_face = fe_face->get_JxW();
                const std::vector<Point> & face_normals = fe_face->get_normals();
                const std::vector<Point> & qface_point = fe_face->get_xyz();
                fe_face->reinit(elem, side);
                
                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    Real g = problem_.computeGradient(stateAdj, qface_point[qp]) + actual_lagrange_;
                    
                    for ( Index k = 0; k < gradJ_.cols(); ++k )
                    {
                        for ( Index l = 0; l < gradJ_.rows(); ++l )
                        {
                            Real b = basisFunction( psi( reference_nodes_(count * quadNodesNo + qp) ), k, l);
                            
                            for ( Index i = 0; i < mesh_->mesh_dimension(); ++i )
                            {
                                gradJ_(gradJ_.rows() - l - 1, k)(i) += b * g * JxW_face[qp] * (boundingBox_.second(i) - boundingBox_.first(i)) * face_normals[qp](i);
                            }
                        }
                    }
                }
                
                ++count;
            }
        }
    }
}

void FFD::applyPerturbation(const EquationSystems & perturbation)
{
    for ( Index k = 0; k < mu_.cols(); ++k )
    {
        for ( Index l = 0; l < mu_.rows(); ++l )
        {
            for ( Index i = 0; i < mesh_->mesh_dimension(); ++i )
            {
                mu_(l, k)(i) -= step_ * gradJ_(l, k)(i);
            }
        }
    }
    
    // Fix control points.
    problem_.fixCP(CP_grid_, mu_);
    
    Mesh::const_element_iterator       ref_el     = reference_mesh_.active_local_elements_begin();
    const Mesh::const_element_iterator ref_end_el = reference_mesh_.active_local_elements_end();
    
    Mesh::const_element_iterator       el     = mesh_->active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh_->active_local_elements_end();
    
    std::vector<bool> hasMoved(mesh_->n_nodes(), false);
    
    for ( ; ref_el != ref_end_el, el != end_el ; ++ref_el, ++el )
    {
        const Elem * ref_elem = *ref_el;
        const Elem * elem     = *el;
        
        Node * ref_node;
        Node * node;
        
        Index subPerSide = elem->n_nodes() / 3 - 1;
        
        for ( Index n = 0; n < elem->n_vertices(); ++n ) // loop over vertices
        {
            ref_node = ref_elem->get_node(n);
            node     = elem->get_node(n);
            
            if ( !hasMoved[node->id()] && problem_.toBeMoved(*node) )
            {
                (*node) = deform( (*ref_node) );
                
                hasMoved[node->id()] = true;
            }
        }
        
        for ( Index n = elem->n_vertices(); n < elem->n_nodes(); ++n ) // loop over non-vertices
        {
            node = elem->get_node(n);
            
            if ( !hasMoved[node->id()] && problem_.toBeMoved(*node) )
            {
                Index idA = (n - 3) / subPerSide;              // ID del vertice precedente.
                Index idB = ((n - 3 + 1) / subPerSide) % 3;    // ID del vertice successivo.
                
                Point * node1 = elem->get_node(idA);
                Point * node2 = elem->get_node(idB);
                
                for ( Index c = 0; c < mesh_->mesh_dimension(); ++c )
                {
                    (*node)(c) = 1.0 / (subPerSide + 1) * ((n - 3) - subPerSide * idA + 1) * ( (*node1)(c) + (*node2)(c) ); // 1/(#nodi per lato+1) * (che ordinamento ha il nodo nel lato) * (incremento)
                }
                
                hasMoved[node->id()] = true;
            }
        }
    }
    
    //mesh_->write("DeformedMesh.vtu");
}

Real FFD::basisFunction(const Point & point, const Index & k, const Index & l) const
{
    Index K = CP_grid_.cols() - 1;
    Index L = CP_grid_.rows() - 1;
    
    return ( binomial_coefficient<Real>(K, k) * std::pow(1 - point(0), K - k) * std::pow(point(0), k) *
             binomial_coefficient<Real>(L, l) * std::pow(1 - point(1), L - l) * std::pow(point(1), l) );
}

Point FFD::psi(const Point & point) const
{
    Point ref_point;
    
    for ( Index i = 0; i < mesh_->mesh_dimension(); ++i )
    {
        if ( boundingBox_.second(i) - boundingBox_.first(i) )
        {
            ref_point(i) = (point(i) - boundingBox_.first(i)) / (boundingBox_.second(i) - boundingBox_.first(i));
        }
    }
    
    return ref_point;
}

Point FFD::psiInv(const Point & ref_point) const
{
    Point point;
    
    for ( Index i = 0; i < 3; ++i )
    {
        point(i) = (boundingBox_.second(i) - boundingBox_.first(i)) * ref_point(i) + boundingBox_.first(i);
    }
    
    return point;
}

Point FFD::deform(const Point & point) const
{
    Point deformed_point(point);
    
    for ( Index k = 0; k < CP_grid_.cols(); ++k )
    {
        for ( Index l = 0; l < CP_grid_.rows(); ++l )
        {
            for ( Index i = 0; i < mesh_->mesh_dimension(); ++i )
            {
                deformed_point(i) += basisFunction(psi(point), k, l) * (boundingBox_.second(i) - boundingBox_.first(i)) * mu_(mu_.rows() - l - 1, k)(i);
            }
        }
    }
    
    return deformed_point;
}
