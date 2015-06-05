#include "ShapeOptimization.h"

ShapeOptimization::ShapeOptimization(const Problem & problem, const std::string & directory, const Real & step, const Index & maxIterationsNo, const Real & tolerance, const bool & volume_constraint, const Real & armijoSlope)
    : problem_(problem), plotName_(directory + "/" + problem_.get_name()), mesh_(problem_.get_mesh()), step_(step), maxIterationsNo_(maxIterationsNo), tolerance_(tolerance), volume_constraint_(volume_constraint), armijoSlope_(armijoSlope), initialVolume_(getVolume())
{}

void ShapeOptimization::apply()
{
    std::string nameOutput(plotName_ + "_Output.txt");
    std::ofstream f_out(nameOutput);
    
    std::shared_ptr<EquationSystems> stateAdj;
    std::shared_ptr<EquationSystems> perturbation;
    
    Real costFunction = 0.0;
    Real costFunctionOld = 0.0;
    
    // Armijo's rule.
    bool approved = false;
    std::shared_ptr<Mesh> meshOld;
    
    // Save the reference mesh.
    mesh_->write(plotName_ + "_ReferenceMesh.vtu");
    
    std::cout << std::endl << "Initial volume: " << getVolume() << std::endl << std::endl;
    
    //Resolve
    Index i = 1;
    
    for ( ; i <= maxIterationsNo_; ++i )
    {
        std::cout << "********** Iteration: " << i << " **********" << std::endl << std::endl;
        
        std::string stateAdj_name = plotName_ + "_StateAndAdjoint" + std::to_string(i) + ".vtk";
        std::string perturbation_name = plotName_ + "_Perturbation" + std::to_string(i) + ".vtk";
        
        if ( i == 1 || !approved )
        {
            stateAdj = std::shared_ptr<EquationSystems>(new EquationSystems(*mesh_));
            problem_.resolveStateAndAdjointEquation(*stateAdj, i);
            
            VTKIO (*mesh_).write_equation_systems(stateAdj_name, *stateAdj);
            
            costFunctionOld = problem_.evaluateCostFunction(*stateAdj);
        }
        
        meshOld = std::shared_ptr<Mesh>(new Mesh(*mesh_));
        
        if ( volume_constraint_ )
        {
            if ( i == 1 )
            {
                old_lagrange_ = problem_.lagrangeMult(*stateAdj);
            }
            
            updateLagrange(problem_.lagrangeMult(*stateAdj));
            std::cout << "Lagrange multiplier = " << actual_lagrange_ << std::endl << std::endl;
        }
        
        Mesh mesh_perturbation(*mesh_);
        
        if ( problem_.get_name() == "Elasticity" )
        {
            perturbation = std::shared_ptr<EquationSystems>(new EquationSystems(*mesh_));
        }
        else if ( problem_.get_name() == "StokesEnergy" )
        {
            perturbation = std::shared_ptr<EquationSystems>(new EquationSystems(mesh_perturbation));
        }
        
        std::cout << "Computing the identity perturbation" << std::endl;
        computePerturbation(*perturbation, *stateAdj);
        
        VTKIO (mesh_perturbation).write_equation_systems (perturbation_name, *perturbation);
        
        std::cout << "Deforming the mesh" << std::endl;
        applyPerturbation(*perturbation);
        std::cout << "    Done." << std::endl << std::endl;
        
        try
        {
            checkDomain();
        }
        catch ( const std::exception & genericException )
        {
            throw;
        }
        
        std::cout << "Deformed volume: " << getVolume() << std::endl << std::endl;
        
        std::cout << "Printing the results" << std::endl << std::endl;
        
        std::string name = plotName_ + "_Deformed" + std::to_string(i) + ".vtu";
        mesh_->write(name);
        
        // Armijo's rule.
        std::cout << "*** Armijo's rule ***" << std::endl << std::endl;
        Real gradJ2 = problem_.sqrGradient(*stateAdj);
        
        stateAdj = std::shared_ptr<EquationSystems>(new EquationSystems(*mesh_));
        problem_.resolveStateAndAdjointEquation(*stateAdj, i);
        
        costFunction = problem_.evaluateCostFunction(*stateAdj);
        
        std::cout << std::endl << "New cost function = " << costFunction << std::endl << std::endl;
        
        std::cout << "gradJ2 = " << gradJ2 << std::endl;
        std::cout << "Cost_old - c * step * gradJ2 = " << costFunctionOld - armijoSlope_ * step_ * gradJ2 << std::endl << std::endl;
        
        // Se il funzionale costo non è diminuito.
        if ( costFunction > costFunctionOld - armijoSlope_ * step_ * gradJ2 )
        {
            approved = false;
            
            step_ /= 2.0;
            
            std::cout << "Step updated! New step = " << step_ << std::endl << std::endl;
            
            mesh_ = std::move(meshOld);
        }
        else
        {
            approved = true;
            
            std::cout << "Armijo APPROVED!" << std::endl << std::endl;
            
            // Overwrite output file.
            VTKIO (*mesh_).write_equation_systems(stateAdj_name, *stateAdj);
            
            f_out << i << ", " << costFunction << ", " << std::min(1.0, std::abs(costFunction - costFunctionOld) / costFunctionOld) << ";" << std::endl;
            
            // Stopping criterion.
            if ( std::abs(costFunction - costFunctionOld) <= tolerance_ * costFunctionOld )
            {
                std::cout << "Convergence achieved!!!" << std::endl << std::endl;
                
                ++i;
                
                break;
            }
            
            costFunctionOld = costFunction;
        }
    }
    
    f_out.close();
    
    --i;
    
    // Export ParaView Data file for time-series visualization.
    std::ofstream pvd_out(plotName_ + "_TimeSeries_" + std::to_string(i) + ".pvd");
    
    pvd_out << "<?xml version=\"1.0\"?>" << std::endl;
    pvd_out << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
    pvd_out << "    <Collection>" << std::endl;
    pvd_out << "        <DataSet timestep=\"0\" file=\"" << problem_.get_name() << "_ReferenceMesh_0.vtu\"/>" << std::endl;
    
    for ( Index k = 1; k <= i; ++k )
    {
        pvd_out << "        <DataSet timestep=\"" << k << "\" file=\"" << problem_.get_name() << "_Deformed" << k << "_0.vtu\"/>" << std::endl;
    }
    
    pvd_out << "    </Collection>" << std::endl;
    pvd_out << "</VTKFile>" << std::endl;
    
    pvd_out.close();
}

void ShapeOptimization::updateLagrange(const Real & lagrange)
{
    actual_lagrange_ = 0.5 * (old_lagrange_ + lagrange) + (getVolume() - initialVolume_) / initialVolume_;
    old_lagrange_ = actual_lagrange_;
}

Real ShapeOptimization::getVolume() const
{
    Mesh::const_element_iterator       el     = mesh_->active_local_elements_begin();
    const Mesh::const_element_iterator end_el = mesh_->active_local_elements_end();
    
    Real volume = 0.0;
    
    for ( ; el != end_el ; ++el )
    {
        const Elem* elem = *el;
        volume += elem->volume();
    }
    
    return volume;
}

void ShapeOptimization::checkDomain() const
{
    bool isReversed = false;
    
    dof_id_type nelem = mesh_->n_elem();
    
    for (dof_id_type i = 0; i < nelem && !isReversed; ++i)
    {
        if ( mesh_->elem(i)->volume() <= 0.0 )
        {
            isReversed = true;
        }
    }
    
    if ( isReversed )
    {
        throw std::runtime_error("checkDomain(): the deformed mesh has negative volumes.");
    }
}
