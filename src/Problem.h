/* C++ */

/**
 * @file   Problem.h
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

#ifndef PROBLEM_H
#define PROBLEM_H

#include "typedefs.h"

/**
 * @class Problem
 *
 * @brief Classe astratta comune a tutti i problemi su cui si applica l'ottimizzazione
 *
 * */

class Problem
{
    public:
        /**
         * @brief Costruttore
         * @param[in] mesh : puntatore alla mesh sul quale è definito il problema
         *
         */
        Problem(Mesh);
        
        /**
         * @brief Distruttore (defaulted)
         *
         */
        virtual ~Problem() = default;
        
        /**
         * @brief Metodo astratto per risolvere lo stato e l'aggiunto
         * @param[out] stateAdj       : Sistema d'equazioni che conterrà lo stato e l'aggiunto
         * @param[in] maxIterationsNo : Numero massimo di iterazioni
         *
         */
        virtual void resolveStateAndAdjointEquation(EquationSystems & stateAdj, const Index & maxIterationsNo) const = 0;
        
        /**
         * @brief Metodo astratto per calcolare il valore del funzionale costo
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @return  il valore del funzionale costo
         *
         */
        virtual Real evaluateCostFunction(EquationSystems & stateAdj) const = 0;
        
        /**
         * @brief Metodo astratto per calcolare il valore del gradiente del funzionale costo in un punto
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @param[in] p        : Punto in cui calcolare il gradiente
         * @return  il valore del gradiente nel punto
         *
         */
        virtual Real computeGradient(EquationSystems & stateAdj, const Point & p) const = 0;
        
        /**
         * @brief Metodo astratto per calcolare la norma @f$ L^2 @f$ del gradiente
         * @param[in] stateAdj : Sistema d'equazioni che contiene lo stato e l'aggiunto
         * @return il valore della norma @f$ L^2 @f$ del gradiente
         *
         */
        virtual Real sqrGradient(EquationSystems & stateAdj) const = 0;
        
        // Auxiliary methods.
        
        /**
         * @brief Metodo astratto per calcolare l'estensione armonica qualora fosse previsto dalla tecnica di ottimizzazione di forma
         * @param[out] perturbation   : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
         * @param[in] stateAdj        : Sistema d'equazioni contenente stato e aggiunto
         * @param[in] lagrange        : lagrangiano
         *
         */
        virtual void harmonicExtension(EquationSystems & perturbation, EquationSystems & stateAdj, const Real & lagrange) const = 0;
        
        /**
         * @brief Metodo astratto per valutare se un nodo può essere spostato o deve rimanere fisso
         * @param[in] node : nodo sul quale si vuole avere l'informazione
         * @return vero se il nodo può essere spostato
         *
         */
        virtual bool toBeMoved(const Node & node) const = 0;
        
        /**
         * @brief Metodo astratto per vincolare l'eventuale spostamento di control point
         * @param[in] CP_grid : la griglia dei control point
         * @param[in] mu      : matrice contenente lo spostamento di ciascun control point nelle due direzioni
         *
         */
        virtual void fixCP(const MatrixXp & CP_grid, MatrixXp & mu) const = 0;
        
        /**
         * @brief Metodo astratto che calcola il moltiplicatore di lagrange
         * @param[in] stateAdj : Sistema d'equazioni contenente stato e aggiunto
         * @return il moltiplicatore di lagrange, che in particolare è l'integrale
         * del gradiente moltiplicato per lo spostamento in direzione normale sul
         * bordo diviso l'integrale degli spostamenti in direzione normale sul bordo
         *
         */
        virtual Real lagrangeMult(EquationSystems & stateAdj) const = 0;
        
        /**
         * @brief Restituisce un puntatore alla mesh del problema
         * @return il puntatore alla mesh
         *
         */
        inline std::shared_ptr<Mesh> get_mesh() const;
        
        /**
         * @brief Restituisce il nome del sistema
         * @return il nome del sistema
         *
         */
        inline std::string get_name() const;
        
    protected:
        std::shared_ptr<Mesh> mesh_; /**< @brief puntatore alla mesh su cui è definito il problema */
        
        std::string name_;           /**< @brief nome del problema che si vuole risolvere */
};

inline std::shared_ptr<Mesh> Problem::get_mesh() const
{
    return mesh_;
}

inline std::string Problem::get_name() const
{
    return name_;
}

#endif /* PROBLEM_H */
