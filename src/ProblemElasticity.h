/* C++ */

/**
 * @file   ProblemElasticity.h
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

#ifndef PROBLEMELASTICITY_H
#define PROBLEMELASTICITY_H

#include "Problem.h"
/**
 * @class ProblemElasticity
 * @brief Classe che eredita da Problem e che rappresenta il problema dell'elasticità lineare
 *
 */

class ProblemElasticity : public Problem
{
    public:
        friend class ElasticityHE;
        friend class ElasticityState;
        
        
        /**
         * @brief Costruttore
         * @param[in] mesh   : puntatore alla mesh sul quale è definito il problema
         * @param[in] lambda : Coefficiente di Lamé @f$ \lambda @f$
         * @param[in] mu     : Coefficiente di Lamé @f$ \mu @f$
         *
         */
        ProblemElasticity(Mesh, const Real &, const Real &);
        
        /**
         * @brief Metodo per risolvere lo stato e l'aggiunto
         * @param[out] stateAdj : Sistema d'equazioni che conterrà lo stato e l'aggiunto
         * @param[in]  n        : Numero massimo di iterazioni
         *
         */
        virtual void resolveStateAndAdjointEquation(EquationSystems &, const Index &) const;
        
        /** @brief Calcola il funzionale costo @f$ \int_{1} u_y d\sigma @f$
         * @param[out] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @return il valore del funzionale costo
         */
        virtual Real evaluateCostFunction(EquationSystems &) const;
        
        /**
         * @brief Metodo per calcolare il valore del gradiente del funzionale costo in un punto
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @param[in] p        : Punto in cui calcolare il gradiente
         * @return  il valore del gradiente nel punto
         *
         * In questo caso il valore di @f$ - 2 \mu | \epsilon (u) |^2 - \lambda (tr(\epsilon (u))^2) @f$ in un punto
         *
         */
        virtual Real computeGradient(EquationSystems &, const Point &) const;
        
        /**
         * @brief Metodo per calcolare la norma @f$ L^2 @f$ del gradiente
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @return il valore della norma @f$ L^2 @f$ del gradiente
         *
         */
        virtual Real sqrGradient(EquationSystems &) const;
        
        // Auxiliary functions.
        /** @copydoc Problem::harmonicExtension(EquationSystems &, EquationSystems &, const Real &) const; */
        virtual void harmonicExtension(EquationSystems &, EquationSystems &, const Real &) const;
        /** @copydoc Problem::toBeMoved(const Node &) */
        virtual bool toBeMoved(const Node &) const;
        /** @copydoc Problem::fixCP(const MatrixXp &, MatrixXp &) */
        virtual void fixCP(const MatrixXp &, MatrixXp &) const;
        /** @copydoc Problem::lagrangeMult(EquationSystems &) */
        virtual Real lagrangeMult(EquationSystems &) const;
        
    protected:
        Real coeff_lambda_;   /**< @brief Coefficiente di Lamé @f$ \lambda @f$ */
        Real coeff_mu_;       /**< @brief Coefficiente di Lamé @f$ \mu @f$ */
};

/**
 * @class ElasticityHE
 * @brief Classe contenente i metodi necessari per calcolare l'estensione armonica nel problema dell'elasticità
 *
 */

class ElasticityHE : public System::Assembly
{
    public:
        /**
         * @brief Costruttore
         * @param[in] perturbation : Riferimento ai sistemi d'equazioni per lo spostamento della mesh
         * @param[in] stateAdj     : Riferimento ai sistemi d'equazioni contenente lo stato e l'aggiunto
         * @param[in] lagrange     : Moltiplicatore di lagrange
         * @param[in] problem      : Riferimento costante al problema dell'elasticità
         *
         */
        ElasticityHE(EquationSystems &, EquationSystems &, const Real &, const ProblemElasticity &);
        
        void assemble(); /**< @brief Assembla le matrici e i vettori per calcolare l'estensione armonica nel problema dell'elasticità */
        
    private:
        EquationSystems & perturbation_; /**< @brief Sistemi d'equazioni per gestire lo spostamento della mesh */
        EquationSystems & stateAdj_;     /**< @brief Sistemi d'equazioni contenente lo stato e l'aggiunto */
        Real lagrange_;                  /**< @brief Moltiplicatore di lagrange */
        
        const ProblemElasticity & problem_; /**< @brief Riferimento costante al problema dell'elasticità */
};

/**
 * @class ElasticityState
 * @brief Classe contenente i metodi necessari per calcolare lo stato (il sistema è autoaggiunto) nel problema dell'elasticità
 *
 */
class ElasticityState : public System::Assembly
{
    public:
        /**
         * @brief Costruttore
         * @param[in] stateAdj     : Riferimento ai sistemi d'equazioni contenente lo stato e l'aggiunto
         * @param[in] problem      : Riferimento costante al problema dell'elasticità
         *
         */
        ElasticityState(EquationSystems &, const ProblemElasticity &);
        
        void assemble(); /**< @brief Assembla le matrici e i vettori per calcolare lo stato nel problema dell'elasticità */
        
        /**
         * @brief Calcola il valore del tensore elasticità negli indici desiderati
         * @param[in] i : Primo indice del tensore
         * @param[in] j : Secondo indice del tensore
         * @param[in] k : Terzo indice del tensore
         * @param[in] l : Quarto indice del tensore
         * @return @f$ D_{i,j,k,\ell} = \lambda \delta_{ij} \delta_{k \ell} + \mu (\delta_{ik} \delta_{j \ell} + \delta_{i \ell} \delta_{jk})@f$
         *
         */
        Real evaluateElasticityTensor(Index i, Index j, Index k, Index l) const;
        
    private:
        EquationSystems & stateAdj_;    /**< @brief Sistemi d'equazioni contenente lo stato e l'aggiunto */
        
        const ProblemElasticity & problem_; /**< @brief Riferimento costante al problema dell'elasticità */
};

#endif /* PROBLEMELASTICITY_H */
