/* C++ */

/**
 * @file   BoundaryDisplacement.h
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

#ifndef BOUNDARYDISPLACEMENT_H
#define BOUNDARYDISPLACEMENT_H

#include "ShapeOptimization.h"

/** @class BoundaryDisplacement
 * @brief Classe che eredita da ShapeOptimization. Utilizza la tecnica del boundary local displacement per eseguire l'ottimizzazione di forma.
 *
 *
 */

class BoundaryDisplacement : public ShapeOptimization
{
    public:
        /**
          * @brief Costruttore
          * @param[in] problem           : Problema sul quale si vuole applicare la Shape Optimization
          * @param[in] directory         : Directory in cui salvare i file di output
          * @param[in] step              : Passo iniziale per il metodo di discesa del gradiente
          * @param[in] maxIterationsNo   : Numero massimo di iterazioni
          * @param[in] tolerance         : Tolleranza per il test d'arresto dell'incremento relativo
          * @param[in] volume_constraint : Specifica se applicare o meno il vincolo di volume
          * @param[in] armijoSlope       : Coefficiente di rilassamento per la regola di Armijo.
          *
          */
        BoundaryDisplacement(const Problem &, const std::string &, const Real &, const Index &, const Real &, const bool &, const Real & = 1.0e-4);
        
        /**
         * @brief Calcola la deformazione della mesh
         * @param[out] perturbation    : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
         * @param[in]  stateAdj        : Sistema d'equazioni contenente stato e aggiunto
         *
         */
        virtual void computePerturbation(EquationSystems &, EquationSystems &);
        
        /**
         * @brief Applica la deformazione alla mesh
         * @param[in] perturbation : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
         *
         */
        virtual void applyPerturbation(const EquationSystems &);
};

#endif /* BOUNDARYDISPLACEMENT_H */
