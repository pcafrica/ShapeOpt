/* C++11 */

#include "src/ShapeOptimizationBase.h"

std::string full_path(const std::string &, const std::string &);

/**
 * @brief The @b main function.
 */
int main(const int argc, const char * const * argv, const char * const * envp)
{
    try
    {
        GetPot commandLine(argc, (char **) argv);
        
        const std::string config_directory =
            std::string(commandLine.follow("../config", 2, "-d", "--directory")) + "/";
            
        const std::string config_file =
            commandLine.follow("config.pot", 2, "-f", "--file");
            
        GetPot config(full_path(config_file, config_directory).c_str());
        
        /**
         * Read general parameters.
         */
        const std::string mesh_filename =
            full_path(config("mesh", "cantilever.msh"), config_directory);
            
        std::string directory = config("output_directory", "Elasticity");
        
        const std::string problemName = config("problem", "Elasticity");
        
        const std::string techniqueName =
            config("technique", "BoundaryDisplacement");
            
        std::cout << "Problem: " << problemName << std::endl
                  << "Technique: " << techniqueName << std::endl
                  << std::endl;
                  
        /**
         * Read problem-related parameters.
         */
        const Real lambda = config("Problem/Elasticity/lambda", 13.0);
        const Real mu     = config("Problem/Elasticity/mu", 5.5);
        
        const Real ux = config("Problem/StokesEnergy/ux", 4.0);
        const Real uy = config("Problem/StokesEnergy/uy", 0.0);
        
        /**
         * Read technique-related parameters.
         */
        const Real step = config("Technique/step", 0.125);
        const Real maxIterationsNo = config("Technique/maxIterationsNo", 80);
        const Real tolerance = config("Technique/tolerance", 1.0e-3);
        
        const bool volume_constraint =
            config("Technique/volume_constraint", true);
            
        const Real armijoSlope = config("Technique/armijoSlope", 1.0e-2);
        
        if ( config.vector_variable_size("Technique/boundingBoxSW") != 2
                || config.vector_variable_size("Technique/boundingBoxNE") != 2 )
        {
            throw std::runtime_error("ERROR: wrong bounding box"
                                     " set in the configuration file.");
        }
        
        Point SW(config("Technique/boundingBoxSW", 0.0, 0),
                 config("Technique/boundingBoxSW", 0.0, 1));
                 
        Point NE(config("Technique/boundingBoxNE", 5.0, 0),
                 config("Technique/boundingBoxNE", 4.0, 1));
                 
        std::pair<Point, Point> boundingBox(SW, NE);
        
        const Index subdivisionsX = config("Technique/FFD/subdivisionsX", 4);
        const Index subdivisionsY = config("Technique/FFD/subdivisionsY", 4);
        
        std::pair<Index, Index> subdivisions(subdivisionsX, subdivisionsY);
        
        const Real alpha = config("Technique/FFD_LS/alpha", 0.99);
        const Index order = config("Technique/DesignElement/order", 3);
        
        /**
         * Instantiate problem.
         */
        LibMeshInit init(argc, argv);
        
        std::cout << "Importing the geometry..." << std::endl;
        Mesh mesh(init.comm(), 2);
        mesh.read(mesh_filename);
        mesh.all_second_order();
        std::cout << "    Done." << std::endl;
        
        Problem * problem;
        
        if ( problemName == "Elasticity" )
        {
            problem = new ProblemElasticity(mesh, lambda, mu);
        }
        else if ( problemName == "StokesEnergy" )
        {
            problem = new ProblemStokesEnergy(mesh, ux, uy);
        }
        else
        {
            throw std::runtime_error("ERROR: wrong variable \"problem\""
                                     " set in the configuration file.");
        }
        
        directory = "Plot_" + problem->get_name() + "_" + directory;
        
        if ( system( ("rm -rf " + directory +
                      " && mkdir " + directory + " 2> /dev/null").c_str() ) );
                      
        /**
         * Instantiate technique.
         */
        ShapeOptimization * shapeOptimization;
        
        if ( techniqueName == "BoundaryDisplacement" )
        {
            shapeOptimization = new BoundaryDisplacement(*problem, directory, step, maxIterationsNo, tolerance, volume_constraint, armijoSlope);
        }
        else if ( techniqueName == "FFD" )
        {
            shapeOptimization = new FFD(*problem, directory, step, maxIterationsNo, tolerance, volume_constraint, boundingBox, subdivisions, armijoSlope);
        }
        else if ( techniqueName == "FFD_LS" )
        {
            shapeOptimization = new FFD_LS(*problem, directory, step, maxIterationsNo, tolerance, volume_constraint, boundingBox, subdivisions, alpha, armijoSlope);
        }
        else if ( techniqueName == "DesignElement" )
        {
            shapeOptimization = new DesignElement(*problem, directory, step, maxIterationsNo, tolerance, volume_constraint, boundingBox, order, armijoSlope);
        }
        else
        {
            throw std::runtime_error("ERROR: wrong variable \"technique\""
                                     " set in the configuration file.");
        }
        
        /**
         * Apply.
         */
        shapeOptimization->apply();
        
        /*mesh.write("ReferenceMesh.vtu");
        EquationSystems es(mesh);
        shapeOptimization->applyPerturbation(es);
        problem->get_mesh()->write("DeformedMesh.vtu");*/
        
        delete problem;
        delete shapeOptimization;
    }
    catch ( const std::exception & genericException )
    {
        std::cerr << genericException.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

std::string full_path(const std::string & filename,
                      const std::string & relative_directory)
{
    return (filename[0] == '/')
           ? filename
           : (relative_directory + filename);
}
