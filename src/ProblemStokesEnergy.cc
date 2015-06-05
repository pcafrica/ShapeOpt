#include "ProblemStokesEnergy.h"

ProblemStokesEnergy::ProblemStokesEnergy(Mesh mesh, const Real & ux, const Real & uy)
    : Problem(mesh), ux_(ux), uy_(uy)
{
    name_ = "StokesEnergy";
}

void ProblemStokesEnergy::resolveStateAndAdjointEquation(EquationSystems & stateAdj, const Index & n) const
{
    /*
     * Define state problem.
     */
    LinearImplicitSystem & stateSystem = stateAdj.add_system<LinearImplicitSystem> (name_);
    unsigned int u_var = stateSystem.add_variable ("u", SECOND);
    unsigned int v_var = stateSystem.add_variable ("v", SECOND);
    stateSystem.add_variable ("p", FIRST);
    
    StokesEnergyState assState(stateAdj, *this);
    stateAdj.get_system(name_).attach_assemble_object(assState);
    
    // Inlet.
    {
        std::vector<unsigned int> variables(2);
        variables[0] = u_var;
        variables[1] = v_var;
        
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(1);
        
        StokesEnergyBC dirichlet(u_var, v_var, ux_, uy_);
        
        DirichletBoundary stateDirichlet_bc(boundary_ids, variables, &dirichlet);
        stateAdj.get_system(name_).get_dof_map().add_dirichlet_boundary(stateDirichlet_bc);
    }
    
    // Simmetry.
    {
        std::vector<unsigned int> variables(1);
        variables[0] = v_var;
        
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(3);
        
        ZeroFunction<> dirichlet;
        DirichletBoundary stateDirichlet_bc(boundary_ids, variables, &dirichlet);
        stateAdj.get_system(name_).get_dof_map().add_dirichlet_boundary(stateDirichlet_bc);
    }
    
    // No-slip.
    {
        std::vector<unsigned int> variables(2);
        variables[0] = u_var;
        variables[1] = v_var;
        
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(4);
        
        ZeroFunction<> dirichlet;
        DirichletBoundary stateDirichlet_bc(boundary_ids, variables, &dirichlet);
        stateAdj.get_system(name_).get_dof_map().add_dirichlet_boundary(stateDirichlet_bc);
    }
    
    /*
     * Define adjoint problem.
     */
    LinearImplicitSystem & adjointSystem = stateAdj.add_system<LinearImplicitSystem> (name_ + "Adjoint");
    adjointSystem.add_variable ("au", SECOND);
    adjointSystem.add_variable ("av", SECOND);
    adjointSystem.add_variable ("ap", FIRST);
    
    StokesEnergyAdjoint assAdjoint(stateAdj, *this);
    stateAdj.get_system(name_ + "Adjoint").attach_assemble_object(assAdjoint);
    
    // Inlet and no-slip.
    {
        std::vector<unsigned int> variables(2);
        variables[0] = u_var;
        variables[1] = v_var;
        
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(1);
        boundary_ids.insert(4);
        
        ZeroFunction<> dirichlet;
        
        DirichletBoundary stateDirichlet_bc(boundary_ids, variables, &dirichlet);
        stateAdj.get_system(name_ + "Adjoint").get_dof_map().add_dirichlet_boundary(stateDirichlet_bc);
    }
    
    // Simmetry.
    {
        std::vector<unsigned int> variables(1);
        variables[0] = v_var;
        
        std::set<boundary_id_type> boundary_ids;
        boundary_ids.insert(3);
        
        ZeroFunction<> dirichlet;
        DirichletBoundary stateDirichlet_bc(boundary_ids, variables, &dirichlet);
        stateAdj.get_system(name_ + "Adjoint").get_dof_map().add_dirichlet_boundary(stateDirichlet_bc);
    }
    
    // Initialize.
    stateAdj.init();
    
    /*
     * Solve the state problem.
     */
    std::cout << "Solving the State Equation." << std::endl;
    stateAdj.get_system(name_).solve();
    std::cout << "    Done." << std::endl;
    
    /*
     * Solve the adjoint problem.
     */
    std::cout << "Solving the Adjoint Equation." << std::endl;
    stateAdj.get_system(name_ + "Adjoint").solve();
    std::cout << "    Done." << std::endl;
}

Real ProblemStokesEnergy::evaluateCostFunction(EquationSystems & stateAdj) const
{
    const MeshBase & mesh = stateAdj.get_mesh();
    Mesh::const_element_iterator  el     = mesh.active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh.active_local_elements_end();
    
    FEType fe_type(FIRST, LAGRANGE);
    AutoPtr<FEBase> fe (FEBase::build(2, fe_type));
    QGauss qrule (2, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);
    
    const std::vector<std::vector<Real> >&  phi = fe->get_phi();
    const std::vector<Real>& JxW = fe->get_JxW();
    const std::vector<Point> & xyz = fe->get_xyz();
    
    Real sum = 0.0;
    
    for ( ; el != end_el ; ++el)
    {
        const Elem* elem = *el;
        
        fe->reinit(elem);
        
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            Gradient du(stateAdj.get_system(name_).point_gradient(0, xyz[qp]));
            Gradient dv(stateAdj.get_system(name_).point_gradient(1, xyz[qp]));
            
            for (unsigned int i = 0; i < phi.size(); i++)
            {
                sum += 0.5 * (du * du + dv * dv) * JxW[qp];
            }
        }
    }
    
    return sum;
}

Real ProblemStokesEnergy::computeGradient(EquationSystems & stateAdj, const Point & p) const
{
    Gradient du = stateAdj.get_system(name_).point_gradient(0, p);
    Gradient dv = stateAdj.get_system(name_).point_gradient(1, p);
    
    Gradient dau = stateAdj.get_system(name_ + "Adjoint").point_gradient(0, p);
    Gradient dav = stateAdj.get_system(name_ + "Adjoint").point_gradient(1, p);
    
    return (du * dau + dv * dav - 0.5 * (du * du + dv * dv));
}

Real ProblemStokesEnergy::sqrGradient(EquationSystems & stateAdj) const
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

void ProblemStokesEnergy::harmonicExtension(EquationSystems & perturbation, EquationSystems & stateAdj, const Real & lagrange) const
{
    LinearImplicitSystem & system = perturbation.add_system<LinearImplicitSystem> ("Perturbation");
    unsigned int u_var = system.add_variable("u", SECOND, LAGRANGE);
    unsigned int v_var = system.add_variable("v", SECOND, LAGRANGE);
    
    StokesEnergyHE assPerturbation(perturbation, stateAdj, lagrange, *this);
    perturbation.get_system("Perturbation").attach_assemble_object(assPerturbation);
    
    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(1);
    boundary_ids.insert(2);
    boundary_ids.insert(3);
    
    std::vector<unsigned int> variables(2);
    variables[0] = u_var;
    variables[1] = v_var;
    
    ZeroFunction<> zf;
    DirichletBoundary dirichlet_bc(boundary_ids, variables, &zf);
    
    perturbation.get_system("Perturbation").get_dof_map().add_dirichlet_boundary(dirichlet_bc);
    
    perturbation.init();
    perturbation.get_system("Perturbation").solve();
}

bool ProblemStokesEnergy::toBeMoved(const Node & node) const
{
    //return ( node(0) != -0.5 && node(0) != 1.0 && node(1) != -0.3 && node(1) != 0.3 );
    return true;
}

void ProblemStokesEnergy::fixCP(const MatrixXp & CP_grid, MatrixXp & mu) const
{
    for ( Index k = 0; k < mu.cols(); ++k )
    {
        for ( Index l = 0; l < mu.rows(); ++l )
        {
            if
            (
                (k == 0) || (k == mu.cols() - 1) || (l == 0) || (l == mu.rows() - 1)
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

Real ProblemStokesEnergy::lagrangeMult(EquationSystems & stateAdj) const
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
    
    for ( ; el != end_el; ++el)
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
                
                if ( mesh.boundary_info->has_boundary_id (elem, s, 4) ) // NACA.
                {
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
    }
    
    return (num / den);
}

StokesEnergyHE::StokesEnergyHE(EquationSystems & perturbation, EquationSystems & stateAdj, const Real & lagrange, const ProblemStokesEnergy & problem)
    : perturbation_(perturbation), stateAdj_(stateAdj), lagrange_(lagrange), problem_(problem) {}

void StokesEnergyHE::assemble()
{
    const MeshBase & mesh = stateAdj_.get_mesh();
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
        {
            if (elem->neighbor(side) == NULL)
            {
                const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
                const std::vector<Real>& JxW_face = fe_face->get_JxW();
                const std::vector<Point>& face_normals = fe_face->get_normals();
                const std::vector<Point> & qface_point = fe_face->get_xyz();
                fe_face->reinit(elem, side);
                
                if ( mesh.boundary_info->has_boundary_id (elem, side, 4) ) // NACA.
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
        }
        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    }
}

StokesEnergyState::StokesEnergyState(EquationSystems & stateAdj, const ProblemStokesEnergy & problem)
    : stateAdj_(stateAdj), problem_(problem) {}

void StokesEnergyState::assemble()
{
    const MeshBase & mesh = stateAdj_.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem & system = stateAdj_.get_system<LinearImplicitSystem>(problem_.name_);
    
    const unsigned int u_var = system.variable_number ("u");
    const unsigned int v_var = system.variable_number ("v");
    const unsigned int p_var = system.variable_number ("p");
    
    FEType fe_vel_type = system.variable_type(u_var);
    FEType fe_pres_type = system.variable_type(p_var);
    
    AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
    AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
    
    fe_vel->attach_quadrature_rule (&qrule);
    fe_pres->attach_quadrature_rule (&qrule);
    
    const std::vector<Real>& JxW = fe_vel->get_JxW();
    
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
    
    const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
    
    const DofMap & dof_map = system.get_dof_map();
    
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    
    DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
        Kvu(Ke), Kvv(Ke), Kvp(Ke),
        Kpu(Ke), Kpv(Ke), Kpp(Ke);
        
    DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);
    
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_p;
    
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        
        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_p_dofs = dof_indices_p.size();
        
        fe_vel->reinit  (elem);
        fe_pres->reinit (elem);
        
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
        
        Kuu.reposition (u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition (u_var * n_u_dofs, p_var * n_u_dofs, n_u_dofs, n_p_dofs);
        
        Kvu.reposition (v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition (v_var * n_v_dofs, p_var * n_v_dofs, n_v_dofs, n_p_dofs);
        
        Kpu.reposition (p_var * n_u_dofs, u_var * n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var * n_u_dofs, v_var * n_u_dofs, n_p_dofs, n_v_dofs);
        Kpp.reposition (p_var * n_u_dofs, p_var * n_u_dofs, n_p_dofs, n_p_dofs);
        
        Fu.reposition (u_var * n_u_dofs, n_u_dofs);
        Fv.reposition (v_var * n_u_dofs, n_v_dofs);
        Fp.reposition (p_var * n_u_dofs, n_p_dofs);
        
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                    Kuu(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kup(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](0);
                    
                    
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                    Kvv(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kvp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](1);
                    
                    
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                    Kpu(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](0);
                    
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                    Kpv(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](1);
                    
        } // end of the quadrature point qp-loop
        
        dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop
    
    return;
}

StokesEnergyAdjoint::StokesEnergyAdjoint(EquationSystems & stateAdj, const ProblemStokesEnergy & problem)
    : stateAdj_(stateAdj), problem_(problem) {}

void StokesEnergyAdjoint::assemble()
{
    const MeshBase & mesh = stateAdj_.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();
    LinearImplicitSystem & system = stateAdj_.get_system<LinearImplicitSystem>(problem_.name_ + "Adjoint");
    
    const unsigned int u_var = system.variable_number ("au");
    const unsigned int v_var = system.variable_number ("av");
    const unsigned int p_var = system.variable_number ("ap");
    
    FEType fe_vel_type = system.variable_type(u_var);
    FEType fe_pres_type = system.variable_type(p_var);
    
    AutoPtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));
    AutoPtr<FEBase> fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
    
    fe_vel->attach_quadrature_rule (&qrule);
    fe_pres->attach_quadrature_rule (&qrule);
    
    const std::vector<Point>& xyz = fe_vel->get_xyz();
    
    const std::vector<Real>& JxW = fe_vel->get_JxW();
    
    const std::vector<std::vector<RealGradient> >& dphi = fe_vel->get_dphi();
    
    const std::vector<std::vector<Real> >& psi = fe_pres->get_phi();
    
    const DofMap & dof_map = system.get_dof_map();
    
    DenseMatrix<Number> Ke;
    DenseVector<Number> Fe;
    
    DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kup(Ke),
        Kvu(Ke), Kvv(Ke), Kvp(Ke),
        Kpu(Ke), Kpv(Ke), Kpp(Ke);
        
    DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fp(Fe);
    
    std::vector<dof_id_type> dof_indices;
    std::vector<dof_id_type> dof_indices_u;
    std::vector<dof_id_type> dof_indices_v;
    std::vector<dof_id_type> dof_indices_p;
    
    
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    for ( ; el != end_el; ++el)
    {
        const Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        
        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_u_dofs = dof_indices_u.size();
        const unsigned int n_v_dofs = dof_indices_v.size();
        const unsigned int n_p_dofs = dof_indices_p.size();
        
        fe_vel->reinit  (elem);
        fe_pres->reinit (elem);
        
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
        
        Kuu.reposition (u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition (u_var * n_u_dofs, p_var * n_u_dofs, n_u_dofs, n_p_dofs);
        
        Kvu.reposition (v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition (v_var * n_v_dofs, p_var * n_v_dofs, n_v_dofs, n_p_dofs);
        
        Kpu.reposition (p_var * n_u_dofs, u_var * n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var * n_u_dofs, v_var * n_u_dofs, n_p_dofs, n_v_dofs);
        Kpp.reposition (p_var * n_u_dofs, p_var * n_u_dofs, n_p_dofs, n_p_dofs);
        
        Fu.reposition (u_var * n_u_dofs, n_u_dofs);
        Fv.reposition (v_var * n_u_dofs, n_v_dofs);
        Fp.reposition (p_var * n_u_dofs, n_p_dofs);
        
        for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
        {
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                    Kuu(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    
            for (unsigned int i = 0; i < n_u_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kup(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](0);
                    
                    
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                    Kvv(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
                    
            for (unsigned int i = 0; i < n_v_dofs; i++)
                for (unsigned int j = 0; j < n_p_dofs; j++)
                    Kvp(i, j) += -JxW[qp] * psi[j][qp] * dphi[i][qp](1);
                    
                    
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_u_dofs; j++)
                    Kpu(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](0);
                    
            for (unsigned int i = 0; i < n_p_dofs; i++)
                for (unsigned int j = 0; j < n_v_dofs; j++)
                    Kpv(i, j) += -JxW[qp] * psi[i][qp] * dphi[j][qp](1);
                    
                    
            Gradient du(stateAdj_.get_system(problem_.name_).point_gradient(u_var, xyz[qp]));
            Gradient dv(stateAdj_.get_system(problem_.name_).point_gradient(v_var, xyz[qp]));
            
            for (unsigned int i = 0; i < n_u_dofs; i++)
            {
                Fu(i) += JxW[qp] * du * dphi[i][qp];
            }
            
            for (unsigned int i = 0; i < n_v_dofs; i++)
            {
                Fv(i) += JxW[qp] * dv * dphi[i][qp];
            }
        } // end of the quadrature point qp-loop
        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
    } // end of element loop
    
    return;
}
