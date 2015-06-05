/* C++ */

/**
 * @file   FFD.h
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

#ifndef FFD_H
#define FFD_H

#include "ShapeOptimization.h"

/**
 * @class FFD
 *
 * @brief Classe che eredita dalla classe ShapeOptimization, utilizza il metodo della Free Form Deformation utilizzando come funzioni di base le B-Spline
 *
 */
class FFD : public ShapeOptimization
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
         * @param[in] boundingBox       : Punti a nord est e a sud ovest indicanti il range della bounding box
         * @param[in] sub               : Coppia contenente il numero di intervalli in cui suddividere la base e l'altezza della bounding box
         * @param[in] armijoSlope       : Coefficiente di rilassamento per la regola di Armijo.
         *
         */
        FFD(const Problem &, const std::string &, const Real &, const Index &, const Real &, const bool &, const std::pair<Point, Point> &, const std::pair<Index, Index> &, const Real & = 1.0e-4);
        
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
        
        /**
         * @brief calcola la funzione di base k, l per il punto x
         * @param[in] point : punto in cui calcolare la funzione di base
         * @param[in] k     : indice orizzontale della griglia
         * @param[in] l     : indice verticale della griglia
         * @return il valore della funzione
         *
         */
        virtual Real basisFunction(const Point &, const Index &, const Index &) const;
        
        /**
         * @brief mappa la scatola nel quadrato unitario
         * @param[in] point : punto nella scatola da trasformare
         * @return le cordinate del punto nel quadrato unitario
         *
         */
        Point psi(const Point &) const;
        
        /**
         * @brief mappa il quadrato unitario nel rettangolo di partenza
         * @param[in] ref_point : punto nel quadrato di riferimento da trasformare
         * @return le cordinate del punto nel rettangolo di partenza
         *
         */
        Point psiInv(const Point &) const;
        
        /**
         * @brief applica la deformazione al punto
         * @param[in] point  : punto in cui calcolare la deformazione
         * @return punto deformato
         *
         */
        Point deform(const Point &) const;
        
    protected:
        Mesh reference_mesh_;                       /**< @brief mesh di riferimento */
        
        VectorXp reference_nodes_;                  /**< @brief vettore contenente i nodi del bordo nella mesh di riferimento */
        std::pair<Point, Point> boundingBox_;       /**< @brief coppia contenente i punti nord est e sud ovest che definiscono la scatola */
        std::pair<Index, Index> sub_;               /**< @brief coppia di numeri indicanti il numero di suddivisioni in orizzontale e in verticale */
        MatrixXp CP_grid_;                          /**< @brief matrice contenente i control point */
        MatrixXp mu_;                               /**< @brief matrice contenente gli spostamenti desiderati per i control point */
        MatrixXp gradJ_;                            /**< @brief matrice contenente il gradiente in funzione dei control point */
        
        bool firstTime_;                            /**< @brief booleano: vero se è la prima volta che calcola la perturbazione dell'identità */
};

#endif /* FFD_H */
