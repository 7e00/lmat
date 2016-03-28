#ifndef NARUTOACM_MATOPSMUL_H_
#define NARUTOACM_MATOPSMUL_H_

#include <assert.h>
#include <string.h>
#include "scalar.h"
#include "matshell.h"
#include "matrix.h"
#include "binaryexpr.h"
#include "matop.h"

namespace narutoacm
{

using meta::types::prio_type;
using meta::types::inferior_type;
using meta::element_type;
using meta::rows;
using meta::cols;
using meta::helper::max;
using meta::helper::is_same_type;
using meta::operator_num;
using meta::is_matrix_expr;

class MatOpSMul : public MatOp
{
    friend class MatOp;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
    friend class BinaryExprHelper;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;

public:
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2>
    static BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<eT2>> GenBinaryExpr(const MatrixBase<eT1,ExprT1> &expr1, Scalar<eT2> sca2)
    {
        return BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<eT2>>(expr1.Derived(), sca2);
    }
    template <typename eT, int M, int N, typename eT1, typename eT2, typename ExprT2>
    static BinaryExpr<eT,M,N,MatOpSMul,Scalar<eT1>,ExprT2> GenBinaryExpr(Scalar<eT1> sca1, const MatrixBase<eT2,ExprT2> &expr2)
    {
        return BinaryExpr<eT,M,N,MatOpSMul,Scalar<eT1>,ExprT2>(sca1, expr2.Derived());
    }

protected:

    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static int rows(const BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT>> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static int rows(const BinaryExpr<eT,M,N,MatOpSMul,Scalar<SeT>,ExprT1> &expr)
    {
        return expr.expr2_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static int cols(const BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT>> &expr)
    {
        return expr.expr1_.Cols();
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static int cols(const BinaryExpr<eT,M,N,MatOpSMul,Scalar<SeT>,ExprT1> &expr)
    {
        return expr.expr2_.Cols();
    }

    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static eT eval(const BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT>> &expr, int idx)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, SeT>::result_type T;
        return T(expr.expr1_(idx)) * T(expr.expr2_);
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static eT eval(const BinaryExpr<eT,M,N,MatOpSMul,Scalar<SeT>,ExprT1> &expr, int idx)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, SeT>::result_type T;
        return T(expr.expr1_) * T(expr.expr2_(idx));
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static eT eval(const BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT>> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, SeT>::result_type T;
        return T(expr.expr1_(r, c)) * T(expr.expr2_);
    }
    template <typename eT, int M, int N, typename ExprT1, typename SeT>
    static eT eval(const BinaryExpr<eT,M,N,MatOpSMul,Scalar<SeT>,ExprT1> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, SeT>::result_type T;
        return T(expr.expr1_) * T(expr.expr2_(r, c));
    }

    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2>
    static void smul(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, Scalar<eT2> sca2, Cache &cache, int relation1 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));

        typedef typename prio_type<typename prio_type<eT1, eT>::result_type, eT2>::result_type T;

        if (relation1 == 0)
        {
            if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
            {
                MatOp::eval(res, expr1.Derived(), cache, 0);
                for (int i = 0; i < expr1.Elems(); ++i)
                    res(i) *= sca2;
            }
            else
            {
                res.SetShape(expr1.Rows(), expr1.Cols());
                for (int i = 0; i < expr1.Elems(); ++i)
                {
                    res(i) = T(expr1(i)) * T(sca2);
                }
            }
        }
        else if (relation1 == 2)
        {
            auto block = cache.Acquire();
            MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            smul(ms, expr1.Derived(), sca2, cache, 0);
            MatOp::eval(res, ms, cache, 0);
        }
        else
        {
            res.SetShape(expr1.Rows(), expr1.Cols());
            for (int i = 0; i < expr1.Elems(); ++i)
                res(i) = T(expr1(i)) * T(sca2);
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2>
    static void smul(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, Scalar<eT2> sca2, Cache &cache, int relation1 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));

        if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result
            && relation1 == 0)
        {
            MatOp::eval(res, expr1.Derived(), cache, 0);
            smul(res, res, sca2, cache, 1);
        }
        else
        {
            smul<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value>(expr1.Derived(), sca2, cache);
            auto block = cache.GetResult();
            MatrixShell<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((eT *)(block->ptr), block->size/sizeof(eT));
            ms.SetShape(expr1.Rows(), expr1.Cols());
            MatOp::eval(res, ms, cache, 0);
        }
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2>
    static void smul(const SimpleMatrixBase<eT1,ExprT1> &expr1, Scalar<eT2> sca2, Cache &cache)
    {
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        smul(ms, expr1.Derived(), sca2, cache, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2>
    static void smul(const MatrixBase<eT1,ExprT1> &expr1, Scalar<eT2> sca2, Cache &cache)
    {
        typedef typename prio_type<eT,eT1>::result_type T;

        MatOp::eval<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((T *)(block1->ptr), block1->size/sizeof(T));
        ms1.SetShape(expr1.Rows(), expr1.Cols());

        MatrixShell<eT,M,N> ms((eT *)(block1->ptr), block1->size/sizeof(eT));
        smul(ms, ms1, sca2, cache, 1);
        cache.StoreResult(block1);
    }

    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,MatOpSMul,ExprT2,Scalar<SeT>> &expr, Cache &cache, int relation = -1)
    {
        smul(res, expr.expr1_, expr.expr2_, cache, relation);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename SeT>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,MatOpSMul,Scalar<SeT>,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        smul(res, expr.expr2_, expr.expr1_, cache, relation);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT1, typename SeT>
    static void eval(const BinaryExpr<eT1,M1,N1,MatOpSMul,ExprT1,Scalar<SeT>> &expr, Cache &cache)
    {
        smul<eT,M,N>(expr.expr1_, expr.expr2_, cache);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT1, typename SeT>
    static void eval(const BinaryExpr<eT1,M1,N1,MatOpSMul,Scalar<SeT>,ExprT1> &expr, Cache &cache)
    {
        smul<eT,M,N>(expr.expr2_, expr.expr1_, cache);
    }
    
protected:
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    MatOpSMul(EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpSMul,ExprT2,ExprT3> &expr, int relation = -1)
        : MatOp(mat, expr, relation)
    {
    }

};

template <typename expreT1, typename ExprT1, typename expreT2>
inline const BinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value,
                        MatOpSMul, ExprT1, Scalar<expreT2>> operator*(const MatrixBase<expreT1,ExprT1> &expr1, Scalar<expreT2> sca2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value>
                        (expr1.Derived(), sca2);
}
template <typename expreT1, typename expreT2, typename ExprT2>
inline const BinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT2>::value,
                        meta::cols<ExprT2>::value,
                        MatOpSMul, Scalar<expreT1>, ExprT2> operator*(Scalar<expreT1> sca1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT2>::value,
                        meta::cols<ExprT2>::value>
                        (sca1, expr2.Derived());
}
template <typename expreT1, typename ExprT1>
inline const BinaryExpr<typename prio_type<expreT1,int>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value,
                        MatOpSMul, ExprT1, Scalar<int>> operator*(const MatrixBase<expreT1,ExprT1> &expr1, int sca2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,int>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value>
                        (expr1.Derived(), Scalar<int>(sca2));
}
template <typename expreT1, typename ExprT1>
inline const BinaryExpr<typename prio_type<expreT1,int>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value,
                        MatOpSMul, Scalar<int>, ExprT1> operator*(int sca1, const MatrixBase<expreT1,ExprT1> &expr2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,int>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value>
                        (Scalar<int>(sca1), expr2.Derived());
}
template <typename expreT1, typename ExprT1>
inline const BinaryExpr<typename prio_type<expreT1,double>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value,
                        MatOpSMul, ExprT1, Scalar<double>> operator*(const MatrixBase<expreT1,ExprT1> &expr1, double sca2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,double>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value>
                        (expr1.Derived(), Scalar<double>(sca2));
}
template <typename expreT1, typename ExprT1>
inline const BinaryExpr<typename prio_type<expreT1,double>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value,
                        MatOpSMul, Scalar<double>, ExprT1> operator*(double sca1, const MatrixBase<expreT1,ExprT1> &expr2)
{
    return MatOpSMul::GenBinaryExpr<typename prio_type<expreT1,double>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT1>::value>
                        (Scalar<double>(sca1), expr2.Derived());
}

} // end of namespace narutoacm

#endif

