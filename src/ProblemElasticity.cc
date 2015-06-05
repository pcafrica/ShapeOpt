#include "ProblemElasticity.h"

ProblemElasticity::ProblemElasticity(Mesh mesh, const Real & lambda, const Real & mu)
    : Problem(mesh), coeff_lambda_(lambda), coeff_mu_(mu)
{
    name_ = "Elasticity";
}

void ProblemElasticity::resolveStateAndAdjointEquation(EquationSystems & stateAdj, const Index & n) const
{
    LinearImplicitSystem & system = stateAdj.add_system<LinearImplicitSystem> (name_);
    unsigned int u_var = system.add_variable("u", SECOND, LAGRANGE);
    unsigned int v_var = system.add_variable("v", SECOND, LAGRANGE);
    
    ElasticityState assState(stateAdj, *this);
    stateAdj.get_system(name_).attach_assemble_object(assState);
    
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(3);
    boundary_ids.insert(5);
    
    std::vector<unsigned int> variables(2);
    variables[0] = u_var;
    variables[1] = v_var;
    
    ZeroFunction<> zf;
    DirichletBoundary dirichlet_bc(boundary_ids, variables, &zf);
    
    stateAdj.get_system(name_).get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    
    stateAdj.init();
    
    std::cout << "Solving the State Equation." << std::endl;
    stateAdj.get_system(name_).solve();
    std::cout << "    Done." << std::endl;
}

Real ProblemElasticity::evaluateCostFunction(EquationSystems & stateAdj) const
{
    const MeshBase & mesh = *mesh_;
    
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<FEBase> fe (FEBase::build(2, fe_type));
    QGauss qrule (2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(2, fe_type));
    QGauss qface(1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    Real sum = 0.0;
    
    Mesh::const_element_iterator       el     = mesh.active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el ; ++el)
    {
        const Elem* elem = *el;
        
        for (unsigned int s = 0; s < elem->n_sides(); s++)
        {
            if (elem->neighbor(s) == NULL)
            {
                if ( mesh.boundary_info->has_boundary_id (elem, s, 1) )
                {
                    const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                    const std::vector<Real>& JxW_face = fe_face->get_JxW();
                    const std::vector<Point> & qface_point = fe_face->get_xyz();
                    
                    fe_face->reinit(elem, s);
                    
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        Real value_sol = stateAdj.get_system(name_).point_value(1, qface_point[qp]);
                        
                        for (unsigned int i = 0; i < phi_face.size(); i++)
                        {
                            sum += std::abs((-1.0) * value_sol * JxW_face[qp] * phi_face[i][qp]);
                        }
                    }
                }
            }
        }
    }
    
    return sum;
}

Real ProblemElasticity::computeGradient(EquationSystems & stateAdj, const Point & p) const
{
    Gradient du = stateAdj.get_system(name_).point_gradient(0, p);
    Gradient dv = stateAdj.get_system(name_).point_gradient(1, p);
    
    Real sum;
    sum = -2. * coeff_mu_ * (
              du(0) * du(0) + dv(1) * dv(1) +
              (
                  (dv(1) + du(0)) *
                  (dv(1) + du(0))
              ) / 4 +
              (
                  (du(1) + dv(0)) *
                  (du(1) + dv(0))
              ) / 4.) - coeff_lambda_ * (du(0) + dv(1)) * (du(0) + dv(1));
    return sum;
}

Real ProblemElasticity::sqrGradient(EquationSystems & stateAdj) const
{
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<FEBase> fe (FEBase::build(2, fe_type));
    QGauss qrule (2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(2, fe_type));
    QGauss qface(1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    const MeshBase & mesh = stateAdj.get_mesh();
    
    Mesh::const_element_iterator       el     = mesh.active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh.active_local_elements_end();
    
    Real gradJ2 = 0.0;
    
    for ( ; el != end_el ; ++el)
    {
        const Elem* elem = *el;
        
        for (unsigned int s = 0; s < elem->n_sides(); s++)
        {
            if (elem->neighbor(s) == NULL)
            {
                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                const std::vector<Point> & qface_point = fe_face->get_xyz();
                
                fe_face->reinit(elem, s);
                
                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    Real g = - computeGradient(stateAdj, qface_point[qp]);
                    
                    for (unsigned int i = 0; i < phi_face.size(); i++)
                    {
                        gradJ2 += g * g * JxW_face[qp] * phi_face[i][qp];
                    }
                }
            }
        }
    }
    
    return gradJ2;
}

void ProblemElasticity::harmonicExtension(EquationSystems & perturbation, EquationSystems & stateAdj, const Real & lagrange) const
{
    LinearImplicitSystem & system = perturbation.add_system<LinearImplicitSystem> ("Perturbation");
    unsigned int u_var = system.add_variable("u", SECOND, LAGRANGE);
    unsigned int v_var = system.add_variable("v", SECOND, LAGRANGE);
    
    ElasticityHE assPerturbation(perturbation, stateAdj, lagrange, *this);
    perturbation.get_system("Perturbation").attach_assemble_object(assPerturbation);
    
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(1);
    boundary_ids.insert(3);
    boundary_ids.insert(5);
    
    std::vector<unsigned int> variables(2);
    variables[0] = u_var;
    variables[1] = v_var;
    
    ZeroFunction<> zf;
    DirichletBoundary dirichlet_bc(boundary_ids, variables, &zf);
    
    perturbation.get_system("Perturbation").get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    
    perturbation.init();
    perturbation.get_system("Perturbation").solve();
}

bool ProblemElasticity::toBeMoved(const Node &) const
{
    return true;
}

void ProblemElasticity::fixCP(const MatrixXp & CP_grid, MatrixXp & mu) const
{
    for ( Index k = 0; k < mu.cols(); ++k )
    {
        for ( Index l = 0; l < mu.rows(); ++l )
        {
            if
            (
                (k == 0) || (k == mu.cols() - 1)
            )
            {
                for ( Index i = 0; i < mesh_->mesh_dimension(); ++i )
                {
                    mu(l, k)(i) = 0.0;
                }
            }
        }
    }
}

Real ProblemElasticity::lagrangeMult(EquationSystems & stateAdj) const
{
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<FEBase> fe (FEBase::build(2, fe_type));
    QGauss qrule (2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(2, fe_type));
    QGauss qface(1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    const MeshBase & mesh = *mesh_;
    
    Mesh::const_element_iterator       el     = mesh.active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh.active_local_elements_end();
    
    Real num = 0.0;
    Real den = 0.0;
    
    for ( ; el != end_el ; ++el)
    {
        const Elem* elem = *el;
        
        for (unsigned int s = 0; s < elem->n_sides(); s++)
        {
            if (elem->neighbor(s) == NULL)
            {
                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                const std::vector<Point> & qface_point = fe_face->get_xyz();
                
                fe_face->reinit(elem, s);
                
                for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                {
                    Real f = - computeGradient(stateAdj, qface_point[qp]);
                    
                    for (unsigned int i = 0; i < phi_face.size(); i++)
                    {
                        num += f * JxW_face[qp] * phi_face[i][qp];
                        den +=     JxW_face[qp] * phi_face[i][qp];
                    }
                }
            }
        }
    }
    
    return (num / den);
}

ElasticityHE::ElasticityHE(EquationSystems & perturbation, EquationSystems & stateAdj, const Real & lagrange, const ProblemElasticity & problem)
    : perturbation_(perturbation), stateAdj_(stateAdj), lagrange_(lagrange), problem_(problem) {}

void ElasticityHE::assemble()
{
    const MeshBase & mesh = perturbation_.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem & system = perturbation_.get_system<LinearImplicitSystem>("Perturbation");
    
    const unsigned int u_var = system.variable_number ("u");
    const unsigned int v_var = system.variable_number ("v");
    
    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    
    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kvu(Ke), Kvv(Ke);
    DenseSubVector<Number> Fu(Fe), Fv(Fe);
    
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        
        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        
        fe->reinit (elem);
        
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
        
        Kuu.reposition (u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        
        Kvu.reposition (v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        
        Fu.reposition (u_var * n_u_dofs, n_u_dofs);
        Fv.reposition (v_var * n_u_dofs, n_v_dofs);
        
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < n_u_dofs; ++i)
                for (unsigned int j = 0; j < n_u_dofs; ++j)
                    Kuu(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    
            for (unsigned int i = 0; i < n_v_dofs; ++i)
                for (unsigned int j = 0; j < n_v_dofs; ++j)
                    Kvv(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
        }
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
            if (elem->neighbor(side) == NULL)
            {
                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                const std::vector<Point>& face_normals = fe_face->get_normals();
                const std::vector<Point> & qface_point = fe_face->get_xyz();
                fe_face->reinit(elem, side);
                
                if ( mesh.boundary_info->has_boundary_id (elem, side, 0) || mesh.boundary_info->has_boundary_id (elem, side, 2) || mesh.boundary_info->has_boundary_id (elem, side, 4) ) // Apply a traction on the right side
                {
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        Real g = - problem_.computeGradient(stateAdj_, qface_point[qp]) - lagrange_;
                        
                        for (unsigned int i = 0; i < n_u_dofs; i++)
                        {
                            Fu(i) += g * JxW_face[qp] * face_normals[qp](0) * phi_face[i][qp];
                        }
                        
                        for (unsigned int i = 0; i < n_v_dofs; i++)
                        {
                            Fv(i) += g * JxW_face[qp] * face_normals[qp](1) * phi_face[i][qp];
                        }
                    }
                }
            }
            
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    }
}

ElasticityState::ElasticityState(EquationSystems & stateAdj, const ProblemElasticity & problem)
    : stateAdj_(stateAdj), problem_(problem) {}

void ElasticityState::assemble()
{
    const MeshBase & mesh = stateAdj_.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem & system = stateAdj_.get_system<LinearImplicitSystem>(problem_.name_);
    
    const unsigned int u_var = system.variable_number ("u");
    const unsigned int v_var = system.variable_number ("v");
    
    const DofMap& dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(0);
    AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);
    
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
    
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    
    DenseSubMatrix<Number> Kuu(Ke), Kuv(Ke), Kvu(Ke), Kvv(Ke);
    DenseSubVector<Number> Fu(Fe), Fv(Fe);
    
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        
        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        
        fe->reinit (elem);
        
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
        
        Kuu.reposition (u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        
        Kvu.reposition (v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        
        Fu.reposition (u_var * n_u_dofs, n_u_dofs);
        Fv.reposition (v_var * n_u_dofs, n_v_dofs);
        
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 0, C_k = 0;
                    
                    C_j = 0, C_l = 0;
                    Kuu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 0;
                    Kuu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 0, C_l = 1;
                    Kuu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 1;
                    Kuu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }
                
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                {
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 0, C_k = 1;
                    
                    
                    C_j = 0, C_l = 0;
                    Kuv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 0;
                    Kuv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 0, C_l = 1;
                    Kuv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 1;
                    Kuv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }
                
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                {
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 1, C_k = 0;
                    
                    
                    C_j = 0, C_l = 0;
                    Kvu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 0;
                    Kvu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 0, C_l = 1;
                    Kvu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 1;
                    Kvu(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }
                
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                {
                    unsigned int C_i, C_j, C_k, C_l;
                    C_i = 1, C_k = 1;
                    
                    
                    C_j = 0, C_l = 0;
                    Kvv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 0;
                    Kvv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 0, C_l = 1;
                    Kvv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                    
                    C_j = 1, C_l = 1;
                    Kvv(i, j) += JxW[qp] * (evaluateElasticityTensor(C_i, C_j, C_k, C_l) * dphi[i][qp](C_j) * dphi[j][qp](C_l));
                }
        }
        
        
        for (unsigned int side = 0; side < elem->n_sides(); side++)
            if (elem->neighbor(side) == NULL)
            {
                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                //const std::vector<Point>& face_normals = fe_face->get_normals();
                fe_face->reinit(elem, side);
                
                if ( mesh.boundary_info->has_boundary_id (elem, side, 1) ) // Apply a traction on the right side
                {
                    for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i < n_v_dofs; i++)
                        {
                            Fv(i) += (-1.0) * JxW_face[qp] * phi_face[i][qp];
                        }
                    }
                }
            }
            
            
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    }
}

Real ElasticityState::evaluateElasticityTensor(Index i, Index j, Index k, Index l) const
{
    Real delta_ij = (i == j) ? 1.0 : 0.0;
    Real delta_il = (i == l) ? 1.0 : 0.0;
    Real delta_ik = (i == k) ? 1.0 : 0.0;
    Real delta_jl = (j == l) ? 1.0 : 0.0;
    Real delta_jk = (j == k) ? 1.0 : 0.0;
    Real delta_kl = (k == l) ? 1.0 : 0.0;
    
    return problem_.coeff_lambda_ * delta_ij * delta_kl + problem_.coeff_mu_ * (delta_ik * delta_jl + delta_il * delta_jk);
}