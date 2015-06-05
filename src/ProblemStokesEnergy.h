/* C++ */

/**
 * @file   ProblemStokesEnergy.h
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

#ifndef PROBLEMSTOKESENERGY_H
#define PROBLEMSTOKESENERGY_H

#include "Problem.h"
/**
 * @class ProblemStokesEnergy
 * @brief Classe che eredita da Problem e che rappresenta il problema di Stokes
 *
 */
class ProblemStokesEnergy : public Problem
{
    public:
        friend class StokesEnergyHE;
        friend class StokesEnergyState;
        friend class StokesEnergyAdjoint;
        
        /**
         * @brief Costruttore
         * @param[in] mesh   : puntatore alla mesh sul quale è definito il problema
         * @param[in] ux     : Componente lungo l'asse @f$ x @f$ della velocità in ingresso
         * @param[in] uy     : Componente lungo l'asse @f$ y @f$ della velocità in ingresso
         *
         */
        ProblemStokesEnergy(Mesh, const Real &, const Real &);
        
        /**
         * @brief Metodo per risolvere lo stato e l'aggiunto
         * @param[out] stateAdj : Sistema d'equazioni che conterrà lo stato e l'aggiunto
         * @param[in]  n        : Numero massimo di iterazioni
         *
         */
        virtual void resolveStateAndAdjointEquation(EquationSystems &, const Index &) const;
        
        virtual Real evaluateCostFunction(EquationSystems &) const;
        
        /**
         * @brief Metodo per calcolare il valore del gradiente del funzionale costo in un punto
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @param[in] p        : Punto in cui calcolare il gradiente
         * @return  il valore del gradiente nel punto
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
        Real ux_;   /**< @brief Componente lungo l'asse @f$ x @f$ della velocità in ingresso */
        Real uy_;   /**< @brief Componente lungo l'asse @f$ y @f$ della velocità in ingresso */
};


/**
 * @class StokesEnergyHE
 * @brief Classe contenente i metodi necessari per calcolare l'estensione armonica nel problema di Stokes
 *
 */

class StokesEnergyHE : public System::Assembly
{
    public:
        /**
         * @brief Costruttore
         * @param[in] perturbation : Riferimento ai sistemi d'equazioni per lo spostamento della mesh
         * @param[in] stateAdj     : Riferimento ai sistemi d'equazioni contenente lo stato e l'aggiunto
         * @param[in] lagrange     : Moltiplicatore di lagrange
         * @param[in] problem      : Riferimento costante al problema di Stokes
         *
         */
        StokesEnergyHE(EquationSystems &, EquationSystems &, const Real &, const ProblemStokesEnergy &);
        
        void assemble();    /**< @brief Assembla le matrici e i vettori per calcolare l'estensione armonica nel problema dell'elasticità */
        
    private:
        EquationSystems & perturbation_;    /**< @brief Sistemi d'equazioni per gestire lo spostamento della mesh */
        EquationSystems & stateAdj_;        /**< @brief Sistemi d'equazioni contenente lo stato e l'aggiunto */
        Real lagrange_;                     /**< @brief Moltiplicatore di lagrange */
        
        const ProblemStokesEnergy & problem_;   /**< @brief Riferimento costante al problema di Stokes */
};

/**
 * @class StokesEnergyState
 * @brief Classe contenente i metodi necessari per calcolare lo stato nel problema di Stokes
 *
 */
class StokesEnergyState : public System::Assembly
{
    public:
        /**
         * @brief Costruttore
         * @param[in] stateAdj     : Riferimento ai sistemi d'equazioni contenente lo stato e l'aggiunto
         * @param[in] problem      : Riferimento costante al problema di Stokes
         *
         */
        StokesEnergyState(EquationSystems &, const ProblemStokesEnergy &);
        
        void assemble();    /**< @brief Assembla le matrici e i vettori per calcolare lo stato nel problema di Stokes */
        
    private:
        EquationSystems & stateAdj_;    /**< @brief Sistemi d'equazioni contenente lo stato e l'aggiunto */
        
        const ProblemStokesEnergy & problem_;   /**< @brief Riferimento costante al problema di Stokes */
};

/**
 * @class StokesEnergyAdjoint
 * @brief Classe contenente i metodi necessari per calcolare l'aggiunto nel problema di Stokes
 *
 */
class StokesEnergyAdjoint : public System::Assembly
{
    public:
        /**
         * @brief Costruttore
         * @param[in] stateAdj     : Riferimento ai sistemi d'equazioni contenente lo stato e l'aggiunto
         * @param[in] problem      : Riferimento costante al problema di Stokes
         *
         */
        StokesEnergyAdjoint(EquationSystems &, const ProblemStokesEnergy &);
        
        void assemble();    /**< @brief Assembla le matrici e i vettori per calcolare l'aggiunto nel problema di Stokes*/
        
    private:
        EquationSystems & stateAdj_;    /**< @brief Sistemi d'equazioni contenente lo stato e l'aggiunto */
        
        const ProblemStokesEnergy & problem_;   /**< @brief Riferimento costante al problema di Stokes */
};

/**
 * @class StokesEnergyBC
 * @brief Classe contenente i metodi necessari per imporre le condizioni al bordo nel problema di Stokes
 *
 */
class StokesEnergyBC : public FunctionBase<Number>
{
    public:
        /**
         * @brief Costruttore
         * @param[in] u_var  : Indice della variabile che fa riferimento alla componente @f$ x @f$ della velocità
         * @param[in] v_var  : Indice della variabile che fa riferimento alla componente @f$ y @f$ della velocità
         * @param[in] ux     : Componente lungo l'asse @f$ x @f$ della velocità in ingresso
         * @param[in] uy     : Componente lungo l'asse @f$ y @f$ della velocità in ingresso
         *
         */
        StokesEnergyBC (const Index & u_var, const Index & v_var, const Real & ux, const Real & uy)
            : u_var_(u_var), v_var_(v_var), ux_(ux), uy_(uy)
        {
            this->_initialized = true;
        }
        
        /**
         * @brief Operatore() non implementato
         *
         */
        virtual Number operator() (const Point&, const Real = 0)
        {
            libmesh_not_implemented();
        }
        
        /**
         * @brief Operatore() che restituisce il valore delle condizioni al bordo.
         * @param[in]  p      : Punto del bordo in cui valutare le condizioni al bordo
         * @param[in]  time   : Istante temporale considerato
         * @param[out] output : Vettore contenente le componenti delle condizioni al bordo nel punto e all'istante considerato
         *
         */
        virtual void operator() (const Point& p,
                                 const Real time,
                                 DenseVector<Number>& output)
        {
            output.resize(2);
            output.zero();
            
            output(u_var_) = ux_;
            output(v_var_) = uy_;
        }
        
        /**
         * @brief Metodo per clonare l'oggetto.
         * @returns una nuova copia dell'oggetto attuale.
         *
         */
        virtual AutoPtr<FunctionBase<Number> > clone() const
        {
            return AutoPtr<FunctionBase<Number> > (new StokesEnergyBC(u_var_, v_var_, ux_, uy_));
        }
        
    private:
        const Index u_var_;    /**< @brief Indice della variabile che fa riferimento alla componente @f$ x @f$ della velocità */
        const Index v_var_;    /**< @brief Indice della variabile che fa riferimento alla componente @f$ Y @f$ della velocità */
        const Real ux_;    /**< @brief Componente lungo l'asse @f$ x @f$ della velocità in ingresso */
        const Real uy_;    /**< @brief Componente lungo l'asse @f$ y @f$ della velocità in ingresso */
};

#endif /* PROBLEMSTOKESENERGY_H */
