#include "DesignElement.h"

DesignElement::DesignElement(const Problem & problem, const std::string & directory, const Real & step, const Index & maxIterationsNo, const Real & tolerance, const bool & volume_constraint, const std::pair<Point, Point> & boundingBox, const Index & order, const Real & armijoSlope)
    : ShapeOptimization(problem, directory, step, maxIterationsNo, tolerance, volume_constraint, armijoSlope), reference_mesh_(*problem_.get_mesh()), boundingBox_(boundingBox), firstTime_(true)
{
    // Initialize mu_ and gradJ_.
    mu_    = VectorXr::Zero(2 * order);
    gradJ_ = mu_;
    
    const Real n = mu_.size();
    
    P_ = MatrixXr::Identity(n, n);
    
    for ( Index i = 0; i < P_.rows(); ++i )
    {
        for ( Index j = 0; j < P_.cols(); ++j )
        {
            if ( (i == n / 2 - 1 && j <= n / 2 - 1)
                    || (i == n - 1 && j >= n / 2)
                    || (j == n / 2 - 1 && i <= n / 2 - 1)
                    || (j == n - 1 && i >= n / 2)
               )
            {
                P_(i, j) = -1;
            }
        }
    }
    
    P_(n / 2 - 1, n / 2 - 1) = order - 1;
    P_(n - 1, n - 1) = order - 1;
    
    P_ *= 0.5;
}

void DesignElement::computePerturbation(EquationSystems & perturbation, EquationSystems & stateAdj)
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
    
    Real H = boundingBox_.second(1) - boundingBox_.first(1);
    
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
                    Point p = psi( reference_nodes_(count * quadNodesNo + qp) );
                    
                    Real g = problem_.computeGradient(stateAdj, qface_point[qp]) + actual_lagrange_;
                    
                    for ( Index k = 0; k < gradJ_.size() / 2; ++k )
                    {
                        Real x = p(0);
                        
                        for ( Index i = 0; i < k; ++i )
                        {
                            x *= p(0);
                        }
                        
                        gradJ_(k) += g * p(1) * x * JxW_face[qp] * face_normals[qp](1);
                        gradJ_(gradJ_.size() / 2 + k) += g * (H - p(1)) * x * JxW_face[qp] * face_normals[qp](1);
                    }
                }
                
                ++count;
            }
        }
    }
    
    gradJ_ = P_ * gradJ_;
}

void DesignElement::applyPerturbation(const EquationSystems & perturbation)
{
    mu_ -= step_ * gradJ_;
    
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

Point DesignElement::psi(const Point & point) const
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

Point DesignElement::psiInv(const Point & ref_point) const
{
    Point point;
    
    for ( Index i = 0; i < 3; ++i )
    {
        point(i) = (boundingBox_.second(i) - boundingBox_.first(i)) * ref_point(i) + boundingBox_.first(i);
    }
    
    return point;
}

Point DesignElement::deform(const Point & point) const
{
    Point deformed_point(point);
    
    Real H = boundingBox_.second(1) - boundingBox_.first(1);
    
    Point psiPoint = psi(point);
    
    Real x = psiPoint(0);
    Real y = psiPoint(1);
    
    Real fUp   = mu_(mu_.size() / 2 - 1);
    Real fDown = mu_(mu_.size() - 1);
    
    for ( Index k = 1; k < mu_.size() / 2; ++k )
    {
        fUp *= x;
        fUp += mu_(mu_.size() / 2 - 1 - k);
        
        fDown *= x;
        fDown += mu_(mu_.size() - 1 - k);
    }
    
    fUp   *= x;
    fDown *= x;
    
    deformed_point(1) += H * (y * fUp + (1 - y) * fDown);
    
    return deformed_point;
}
