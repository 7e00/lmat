#ifndef NARUTOACM_MATOPMUL_H_
#define NARUTOACM_MATOPMUL_H_

#include <assert.h>
#include <string.h>
#include "matshell.h"
#include "matrix.h"
#include "binaryexpr.h"
#include "matop.h"

namespace narutoacm
{

using meta::types::prio_type;
using meta::element_type;
using meta::rows;
using meta::cols;
using meta::helper::max;
using meta::helper::is_same_type;

class MatOpMul : public MatOp
{
    friend class MatOp;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
    friend class BinaryExprHelper;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;

public:
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2> GenBinaryExpr(const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2)
    {
        return BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2>(expr1.Derived(), expr2.Derived());
    }

protected:

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int rows(const BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int cols(const BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2> &expr)
    {
        return expr.expr2_.Cols();
    }

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2> &expr, int idx)
    {
#if MATRIXORDER == ROWORDER
        return eval(expr, idx/expr.expr2_.Cols(), idx%expr.expr2_.Cols());
#else
        return eval(expr, idx%expr.expr1_.Rows(), idx/expr.expr1_.Rows());
#endif
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        T res = 0;
        for (int k = 0; k < expr.expr1_.Cols(); ++k)
        {
            res += T(expr.expr1_(r,k)) * T(expr.expr2_(k,c));
        }
        return res;
    }

    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));
        
        typedef typename prio_type<typename prio_type<eT1, eT>::result_type, eT2>::result_type T;

        if (relation1 == 0 && relation2 == 0)
        {
            res.SetShape(expr1.Rows(), expr2.Cols());
            int idx = 0;
#if MATRIXORDER == ROWORDER
            for (int i = 0; i < expr1.Rows(); ++i)
            {
                for (int j = 0; j < expr2.Cols(); ++j)
                {
#else
            for (int j = 0; j < expr2.Cols(); ++j)
            {
                for (int i = 0; i < expr1.Rows(); ++i)
                {
#endif
                    T sum = 0;
                    for (int k = 0; k < expr1.Cols(); ++k)
                    {
                        sum += T(expr1(i,k)) * T(expr2(k,j));
                    }
                    res(idx++) = sum;
                }
            }
        }
        else
        {
            auto block = cache.Acquire();
            MatrixShell<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((eT *)(block->ptr), block->size/sizeof(eT));
            mul(ms, expr1.Derived(), expr2.Derived(), cache, 0, 0);
            MatOp::eval(res, ms, cache, 0);
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));
        
        MatOp::eval<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((eT1 *)(block1->ptr), block1->size/sizeof(eT1));
        ms1.SetShape(expr1.Rows(), expr1.Cols());
        mul(res, ms1, expr2.Derived(), 0, relation2);
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        
        MatOp::eval<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((eT2 *)(block2->ptr), block2->size/sizeof(eT2));
        ms2.SetShape(expr2.Rows(), expr2.Cols());
        mul(res, expr1.Derived(), ms2, relation1, 0);
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        MatOp::eval<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((eT1 *)(block1->ptr), block1->size/sizeof(eT1));
        ms1.SetShape(expr1.Rows(), expr1.Cols());
        MatOp::eval<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((eT2 *)(block2->ptr), block2->size/sizeof(eT2));
        ms2.SetShape(expr2.Rows(), expr2.Cols());
        mul(res, ms1, ms2, cache, 0, 0);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(const SimpleMatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        mul(ms, expr1.Derived(), expr2.Derived(), cache);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(const MatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        MatOp::eval<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((eT1 *)(block1->ptr), block1->size/sizeof(eT1));
        ms1.SetShape(expr1.Rows(), expr1.Cols());
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        mul(ms, ms1, expr2.Derived(), cache, 0, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(const SimpleMatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        MatOp::eval<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((eT2 *)(block2->ptr), block2->size/sizeof(eT2));
        ms2.SetShape(expr2.Rows(), expr2.Cols());
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        mul(ms, expr1.Derived(), ms2, cache, 0, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void mul(const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        MatOp::eval<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((eT1 *)(block1->ptr), block1->size/sizeof(eT1));
        ms1.SetShape(expr1.Rows(), expr1.Cols());
        MatOp::eval<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((eT2 *)(block2->ptr), block2->size/sizeof(eT2));
        ms2.SetShape(expr2.Rows(), expr2.Cols());
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        mul(ms, ms1, ms2, cache, 0, 0);
        cache.StoreResult(block);
    }

    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,MatOpMul,ExprT2,ExprT3> &expr, Cache &cache, int relation = -1)
    {
        relation = (relation == 0 ? 0 : -1);
        mul(res, expr.expr1_, expr.expr2_, cache, relation, relation);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT1, typename ExprT2>
    static void eval(const BinaryExpr<eT1,M1,N1,MatOpMul,ExprT1,ExprT2> &expr, Cache &cache)
    {
        mul<eT,M,N>(expr.expr1_, expr.expr2_, cache);
    }

    
protected:
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    MatOpMul(EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpMul,ExprT2,ExprT3> &expr)
        : MatOp(mat, expr)
    {
    }

};

template <typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
inline const BinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT2>::value,
                        MatOpMul, ExprT1, ExprT2> operator*(const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    assert(expr1.Cols() == expr2.Rows());
    return MatOpMul::GenBinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        meta::rows<ExprT1>::value,
                        meta::cols<ExprT2>::value>
                        (expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm

#endif

