/* C++ */

/**
 * @file   ShapeOptimization.h
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

#ifndef SHAPEOPTIMIZATION_H
#define SHAPEOPTIMIZATION_H

#include "typedefs.h"

#include "Problem.h"
#include "ProblemElasticity.h"
#include "ProblemStokesEnergy.h"

/**
 * @class ShapeOptimization
 *
 * @brief Classe astratta comune a tutte le tecniche di ottimizzazione
 *
 *
 */

class ShapeOptimization
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
        ShapeOptimization(const Problem &, const std::string &, const Real &, const Index &, const Real &, const bool &, const Real & = 1.0e-4);
        
        /**
         * @brief Distruttore (defaulted)
         *
         */
        virtual ~ShapeOptimization() = default;
        
        /**
         * @brief Esegue il ciclo d'ottimizzazione
         *
         */
        void apply();
        
        /**
         * @brief Metodo astratto per calcolare la deformazione della mesh
         * @param[out] perturbation    : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
         * @param[in]  stateAdj : Sistema d'equazioni contenente stato e aggiunto
         *
         */
        virtual void computePerturbation(EquationSystems & perturbation, EquationSystems & stateAdj) = 0;
        
        /**
         * @brief Metodo astratto per applicare la deformazione alla mesh
         * @param[in] perturbation : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
         *
         */
        virtual void applyPerturbation(const EquationSystems & perturbation) = 0;
        
        /**
         * @brief Aggiorna il valore del moltiplicatore di lagrange
         * @f$ l_{k+1} = \frac{ l + l_k}{2} + \frac{ V - V_0 }{ V_0 } @f$
         * @param[in] lagrange : media del gradiente sul bordo
         * @f$ l = \frac{\int_{\partial \Omega}{- \nabla J d \sigma}}{\int_{\partial \Omega}{d \sigma}} @f$
         *
         */
        void updateLagrange(const Real &);
        
        /**
         * @brief Misura l'area della mesh
         * @return il valore dell'area della mesh
         *
         */
        Real getVolume() const;
        
        /**
         * @brief Controlla se non si sono invertiti dei triangoli della mesh in seguito alla deformazione
         *
         */
        void checkDomain() const;
        
    protected:
        const Problem & problem_;       /**< @brief problema che si vuole ottimizzare */
        std::string plotName_;          /**< @brief nome utilizzato nella generazione dei file di output */
        
        std::shared_ptr<Mesh> mesh_;    /**< @brief puntatore alla mesh su cui è definito il problema */
        
        Real step_;                     /**< @brief passo utilizzato per il metodo di discesa del gradiente */
        Index maxIterationsNo_;         /**< @brief numero massimo di iterazioni */
        Real tolerance_;                /**< @brief tolleranza per il test d'arresto dell'incremento relativo */
        bool volume_constraint_;        /**< @brief specifica se applicare o meno il vincolo di volume */
        Real armijoSlope_;              /**< @brief coefficiente di rilassamento per la regola di Armijo */
        
        Real old_lagrange_;             /**< @brief valore del lagrangiano al passo d'ottimizzazione precedente */
        Real actual_lagrange_;          /**< @brief valore del lagrangiano al passo d'ottimizzazione attuale */
        
        Real initialVolume_;            /**< @brief area iniziale della mesh */
};

#endif /* SHAPEOPTIMIZATION_H */
