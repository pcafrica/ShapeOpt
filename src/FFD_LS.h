/* C++ */

/**
 * @file   FFD_LS.h
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

#ifndef FFD_LS_H
#define FFD_LS_H

#include "FFD.h"
/**
 * @class FFD_LS
 *
 * @brief Classe che eredita dalla classe FFD, utilizza il metodo dei minimi quadrati con rilassamento per calcolare gli spostamenti da applicare ai control point
 *
 * La perturbazione dell'identità è rappresentata nel caso bidimensionale da
 * @f$ \vec{\theta}_{FFD} = (\theta_{FFD,x}, \theta_{FFD,y})^T \in \mathbb{R}^2 @f$
 * @f[
 *    \vec{\theta}_{FFD} (\vec{x}, \vec{\mu}) = \sum_{k = 0}^{K} \sum_{\ell = 0}^{L}
 *       b_{k,\ell}^{K, L} ( \vec{\psi} (\vec{x})) \mathfrak{B} \vec{\mu}_{k, \ell}
 * @f]
 * Ora se uniamo la tripla sommatoria e trasformiamo in un vettore @f$ \vec{\mu}_{k, \ell} @f$ e se consideriamo solo
 * i NB nodi del bordo possiamo esprimere, definendo LL = K x L,
 * la precedente formula con delle matrici @f$ B_x, B_y \in \mathbb{R}^{NB \times LL} @f$ ottenendo
 * @f[
 *    \vec{\theta}_{FFD,i} = B_i  \vec{\mu}_i
 * @f]
 */
class FFD_LS : public FFD
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
        * @param[in] alpha             : Parametro di rilassamento per il metodo dei minimi quadrati
        * @param[in] armijoSlope       : Coefficiente di rilassamento per la regola di Armijo.
        *
        */
        FFD_LS(const Problem &, const std::string &, const Real &, const Index &, const Real &, const bool &, const std::pair<Point, Point> &, const std::pair<Index, Index> &, const Real &, const Real & = 1.0e-4);
        
        /**
        * @brief Calcola la deformazione della mesh
        * @param[out] perturbation    : Sistema d'equazioni contenente gli spostamenti da applicare alla mesh
        * @param[in]  stateAdj        : Sistema d'equazioni contenente stato e aggiunto
        *
        */
        virtual void computePerturbation(EquationSystems &, EquationSystems &);
        
    protected:
        Real alpha_;    /**< @brief Parametro di rilassamento per il metodo dei minimi quadrati */
        
        std::vector<Point> border_ref_;   /**< @brief Vettore contenente i punti del bordo*/
        MatrixXr B_x;   /**< @brief Matrice con righe pari al numero di punti sul bordo e colonne pari al numero totale di control point. Rappresenta @f$ \left( b_{k,\ell}^{K, L} ( \vec{\psi} (\vec{x})) \mathfrak{B} \right)_x @f$ per i nodi del bordo */
        MatrixXr B_y;   /**< @brief Matrice con righe pari al numero di punti sul bordo e colonne pari al numero totale di control point. Rappresenta @f$ \left( b_{k,\ell}^{K, L} ( \vec{\psi} (\vec{x})) \mathfrak{B} \right)_y @f$ per i nodi del bordo */
        
        LDLT<MatrixXr> solver_x; /**< @brief Robust Cholesky Decomposition con pivoting della matrice @f$ \alpha B\_x^T B\_x + (1 - \alpha ) I @f$ */
        LDLT<MatrixXr> solver_y; /**< @brief Robust Cholesky Decomposition con pivoting della matrice @f$ \alpha B\_y^T B\_y + (1 - \alpha ) I @f$ */
};

#endif /* FFD_LS_H */
