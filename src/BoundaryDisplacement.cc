#include "BoundaryDisplacement.h"

BoundaryDisplacement::BoundaryDisplacement(const Problem & problem, const std::string & directory, const Real & step, const Index & maxIterationsNo, const Real & tolerance, const bool & volume_constraint, const Real & armijoSlope)
    : ShapeOptimization(problem, directory, step, maxIterationsNo, tolerance, volume_constraint, armijoSlope) {}

void BoundaryDisplacement::computePerturbation(EquationSystems & perturbation, EquationSystems & stateAdj)
{
    problem_.harmonicExtension(perturbation, stateAdj, actual_lagrange_);
}

void BoundaryDisplacement::applyPerturbation(const EquationSystems & perturbation)
{
    Mesh::const_element_iterator       el     = mesh_->active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh_->active_local_elements_end();
    
    std::vector<bool> hasMoved(mesh_->n_nodes(), false);
    
    for ( ; el != end_el ; ++el )
    {
        const Elem * elem = *el;
        Node * node;
        Index subPerSide = elem->n_nodes() / 3 - 1;
        
        for ( Index n = 0; n < elem->n_vertices(); ++n ) // loop over vertices
        {
            node = elem->get_node(n);
            
            if ( !hasMoved[node->id()] && problem_.toBeMoved(*node) )
            {
                for ( Index c = 0; c < mesh_->mesh_dimension(); ++c )
                {
                    (*node)(c) += step_ * (perturbation.get_system("Perturbation").point_value(c, *node));
                }
                
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
}
