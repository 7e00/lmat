#ifndef NARUTOACM_MATOPEMUL_H_
#define NARUTOACM_MATOPEMUL_H_

#include <assert.h>
#include <string.h>
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

class MatOpEMul : public MatOp
{
    friend class MatOp;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
    friend class BinaryExprHelper;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;

public:
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2> GenBinaryExpr(const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2)
    {
        return BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2>(expr1.Derived(), expr2.Derived());
    }

protected:

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int rows(const BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int cols(const BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Cols();
    }

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2> &expr, int idx)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(idx)) * T(expr.expr2_(idx));
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,MatOpEMul,ExprT1,ExprT2> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(r, c)) * T(expr.expr2_(r, c));
    }

    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));

        typedef typename prio_type<typename prio_type<eT1, eT>::result_type, eT2>::result_type T;

        if (relation1 == 0 && relation2 == 0)
        {
            if (is_same_type<eT,T>::result)
            {
                if (operator_num<MatOp,ExprT1>::value >= operator_num<MatOp,ExprT2>::value)
                {
                    MatOp::eval(res, expr1.Derived(), cache, 0);
                    for (int i = 0; i < expr2.Elems(); ++i)
                        res(i) *= expr2(i);
                }
                else
                {
                    MatOp::eval(res, expr2.Derived(), cache, 0);
                    for (int i = 0; i < expr1.Elems(); ++i)
                        res(i) *= expr1(i);
                }
            }
            else if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
            {
                MatOp::eval(res, expr1.Derived(), cache, 0);
                for (int i = 0; i < expr2.Elems(); ++i)
                    res(i) *= expr2(i);
            }
            else if (is_same_type<eT,typename prio_type<eT,eT2>::result_type>::result)
            {
                MatOp::eval(res, expr2.Derived(), cache, 0);
                for (int i = 0; i < expr1.Elems(); ++i)
                    res(i) *= expr1(i);
            }
            else
            {
                typedef typename prio_type<eT1,eT2>::result_type T;
                res.SetShape(expr1.Rows(), expr1.Cols());
                for (int i = 0; i < expr1.Elems(); ++i)
                {
                    res(i) = T(expr1(i)) * T(expr2(i));
                }
            }
        }
        else if (relation1 == 2 || relation2 == 2)
        {
            auto block = cache.Acquire();
            MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            emul(ms, expr1.Derived(), expr2.Derived(), cache, 0, 0);
            MatOp::eval(res, ms, cache, 0);
        }
        else
        {
            res.SetShape(expr1.Rows(), expr1.Cols());
            for (int i = 0; i < expr1.Elems(); ++i)
                res(i) = T(expr1(i)) * T(expr2(i));
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));

        if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result
            && relation1 == 0 && relation2 == 0)
        {
            MatOp::eval(res, expr1.Derived(), cache, 0);
            emul(res, res, expr2.Derived(), cache, 1, 0);
        }
        else
        {
            emul<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value>(expr1.Derived(), expr2.Derived(), cache);
            auto block = cache.GetResult();
            MatrixShell<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((eT *)(block->ptr), block->size/sizeof(eT));
            ms.SetShape(expr1.Rows(), expr1.Cols());
            MatOp::eval(res, ms, cache, 0);
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        emul(res, expr2.Derived(), expr1.Derived(), cache, relation2, relation1);
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));
        
        typedef typename prio_type<typename prio_type<eT1, eT>::result_type,eT2>::result_type T;

        if (relation1 == 0 && relation2 == 0
            && is_same_type<eT,typename prio_type<eT,typename inferior_type<eT1,eT2>::result_type>::result_type>::result)
        {
            if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
            {
                MatOp::eval(res, expr1.Derived(), cache, 0);
                MatOp::eval<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
                auto block2 = cache.GetResult();
                MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((eT2 *)(block2->ptr), block2->size/sizeof(eT2));
                emul(res, res, ms2, cache, 1, 0);
            }
            else// if (is_same_type<eT,typename prio_type<eT,eT2>::result_type>::result)
            {
                MatOp::eval(res, expr2.Derived(), cache, 0);
                MatOp::eval<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
                auto block1 = cache.GetResult();
                MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((eT1 *)(block1->ptr), block1->size/sizeof(eT1));
                emul(res, ms1, res, cache, 0, 1);
            }
        }
        else
        {
            emul<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value>(expr1.Derived(), expr2.Derived(), cache);
            auto block = cache.GetResult();
            MatrixShell<eT,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((eT *)(block->ptr), block->size/sizeof(eT));
            ms.SetShape(expr1.Rows(), expr1.Cols());
            MatOp::eval(res, ms, cache, 0);
        }
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(const SimpleMatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        auto block = cache.Acquire();
        MatrixShell<eT,M,N> ms((eT *)(block->ptr), block->size/sizeof(eT));
        emul(ms, expr1.Derived(), expr2.Derived(), cache, 0, 0);
        cache.StoreResult(block);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(const MatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        typedef typename prio_type<eT,eT1>::result_type T;

        MatOp::eval<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((T *)(block1->ptr), block1->size/sizeof(T));
        ms1.SetShape(expr1.Rows(), expr1.Cols());

        MatrixShell<eT,M,N> ms((eT *)(block1->ptr), block1->size/sizeof(eT));
        emul(ms, ms1, expr2.Derived(), cache, 1, 0);
        cache.StoreResult(block1);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(const SimpleMatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        typedef typename prio_type<eT,eT2>::result_type T;

        MatOp::eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((T *)(block2->ptr), block2->size/sizeof(T));
        ms2.SetShape(expr2.Rows(), expr2.Cols());

        MatrixShell<eT,M,N> ms((eT *)(block2->ptr), block2->size/sizeof(eT));
        emul(ms, expr1.Derived(), ms2, cache, 0, 1);
        cache.StoreResult(block2);
    }
    template <typename eT, int M, int N, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void emul(const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache)
    {
        typedef typename prio_type<eT,typename prio_type<eT1,eT2>::result_type>::result_type T;

        MatOp::eval<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value>(expr1.Derived(), cache);
        auto block1 = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms1((T *)(block1->ptr), block1->size/sizeof(T));
        ms1.SetShape(expr1.Rows(), expr1.Cols());

        MatOp::eval<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value>(expr2.Derived(), cache);
        auto block2 = cache.GetResult();
        MatrixShell<T,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms2((T *)(block2->ptr), block2->size/sizeof(T));
        ms2.SetShape(expr2.Rows(), expr2.Cols());

        MatrixShell<eT,M,N> ms((eT *)(block1->ptr), block1->size/sizeof(eT));
        emul(ms, ms1, ms2, cache, 1, 0);
        cache.StoreResult(block1);
    }

    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,MatOpEMul,ExprT2,ExprT3> &expr, Cache &cache, int relation = -1)
    {
        relation = (relation == 0 ? 0 : -1);
        emul(res, expr.expr1_, expr.expr2_, cache, relation, relation);
    }
    template <typename eT, int M, int N, typename eT1, int M1, int N1, typename ExprT1, typename ExprT2>
    static void eval(const BinaryExpr<eT1,M1,N1,MatOpEMul,ExprT1,ExprT2> &expr, Cache &cache)
    {
        emul<eT,M,N>(expr.expr1_, expr.expr2_, cache);
    }
    
protected:
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    MatOpEMul(EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,MatOpEMul,ExprT2,ExprT3> &expr, int relation = -1)
        : MatOp(mat, expr, relation)
    {
    }

};

template <typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
inline const BinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result,
                        MatOpEMul, ExprT1, ExprT2> ElementMulti(const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    assert(expr1.Rows() == expr2.Rows() && expr1.Cols() == expr2.Cols());
    return MatOpEMul::GenBinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result>
                        (expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm

#endif

