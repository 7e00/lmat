#ifndef NARUTOACM_OPADD_H_
#define NARUTOACM_OPADD_H_

#include <string.h>
#include "matrix.h"
#include "binaryexpr.h"
#include "op.h"

namespace narutoacm
{

using meta::types::prio_type;
using meta::element_type;

class OpAdd : public Op
{
    friend class Op;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExprBase;
    template <typename eT, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;
    template <typename ExprT1, typename ExprT2>
    friend const BinaryExpr<typename prio_type<typename element_type<ExprT1>::type,
                            typename element_type<ExprT2>::type>::result_type,
                            OpAdd, ExprT1, ExprT2> operator+(const MatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2);

protected:
    template <typename ExprT1, typename ExprT2>
    static const BinaryExpr<typename prio_type<typename element_type<ExprT1>::type,
                            typename element_type<ExprT2>::type>::result_type,
                            OpAdd, ExprT1, ExprT2> add(const MatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2)
    {
        return BinaryExpr<typename prio_type<typename element_type<ExprT1>::type,typename element_type<ExprT2>::type>::result_type, OpAdd, ExprT1, ExprT2>(expr1.Derived(), expr2.Derived());
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static int rows(const BinaryExpr<eT,OpAdd,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, typename ExprT1, typename ExprT2>
    static int cols(const BinaryExpr<eT,OpAdd,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Cols();
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,OpAdd,ExprT1,ExprT2> &expr, int idx)
    {
        typedef typename prio_type<typename element_type<ExprT1>::type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(idx)) + T(expr.expr2_(idx));
    }
    template <typename eT, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,OpAdd,ExprT1,ExprT2> &expr, int r, int c)
    {
        typedef typename prio_type<typename element_type<ExprT1>::type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(r, c)) + T(expr.expr2_(r, c));
    }

    template <typename eT, int M1, int N1, int M2, int N2>
    static void add(eT *res, int len, const Matrix<eT,M1,N1> &expr1, const Matrix<eT,M2,N2> &expr2)
    {
        memcpy(res, expr1.Data(), expr1.Elems() * sizeof(eT));
        auto ptr = expr2.Data();
        for (int i = 0; i < expr2.Elems(); ++i)
            res[i] += ptr[i];
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static void add(eT *res, int len, const SimpleMatrixBase<ExprT1> &expr1, const SimpleMatrixBase<ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        Op::eval(res, len, mat1, cache);
        for (int i = 0; i < mat2.Elems(); ++i)
            res[i] += mat2(i);
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static void add(eT *res, int len, const MatrixBase<ExprT1> &expr1, const SimpleMatrixBase<ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        Op::eval(res, len, mat1, cache);
        for (int i = 0; i < mat2.Elems(); ++i)
            res[i] += mat2(i);
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static void add(eT *res, int len, const SimpleMatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        Op::eval(res, len, mat2, cache);
        for (int i = 0; i < mat1.Elems(); ++i)
            res[i] += mat1(i);
    }

    template <typename eT, typename ExprT1, typename ExprT2>
    static void add(eT *res, int len, const MatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        auto block = cache.Acquire();
        Op::eval((eT *)(block->ptr), block->size / sizeof(eT), mat1, cache);
        memcpy(res, block->ptr, mat1.Elems() * sizeof(eT));
        Op::eval((eT *)(block->ptr), block->size / sizeof(eT), mat2, cache);

        eT *ptr = (eT *)(block->ptr);
        for (int i = 0; i < mat2.Elems(); ++i)
            res[i] += ptr[i];
    }

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static void add(Matrix<eT,M,N> &res, const MatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();

        int ar = relation_to_mat(res, mat1);
        int br = relation_to_mat(res, mat2);

        if (ar == 0 && br == 0)
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            add(res.Data(), res.Capacity(), mat1, mat2, cache);
        }
        else if (ar == 0)
        {
            auto block = cache.Acquire();
            Op::eval((eT *)(block->ptr), block->size / sizeof(eT), mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            Op::eval(res.Data(), res.Capacity(), mat1, cache);
            eT *ptr = (eT *)(block->ptr);
            for (int i = 0; i < res.Elems(); ++i)
                res(i) += ptr[i];
        }
        else if (br == 0)
        {
            auto block = cache.Acquire();
            Op::eval((eT *)(block->ptr), block->size / sizeof(eT), mat1, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            Op::eval(res.Data(), res.Capacity(), mat2, cache);
            eT *ptr = (eT *)(block->ptr);
            for (int i = 0; i < res.Elems(); ++i)
                res(i) += ptr[i];
        }
        else
        {
            auto block = cache.Acquire();
            add((eT *)(block->ptr), block->size / sizeof(eT), mat1, mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            memcpy(res.Data(), block->ptr, res.Elems() * sizeof(eT));
        }
    }

    template <typename eT, typename eeT, typename ExprT1, typename ExprT2>
    static void eval(eT *res, int len, const BinaryExpr<eeT,OpAdd,ExprT1,ExprT2> &expr, Cache &cache)
    {
        add(res, len, expr.expr1_, expr.expr2_, cache);
    }
    template <typename eT, int M, int N, typename eeT, typename ExprT1, typename ExprT2>
    static void eval(Matrix<eT,M,N> &res, const BinaryExpr<eeT,OpAdd,ExprT1,ExprT2> &expr, Cache &cache)
    {
        add(res, expr.expr1_, expr.expr2_, cache);
    }
    
protected:
    template <typename eT, int M, int N, typename eeT, typename ExprT1, typename ExprT2>
    OpAdd(Matrix<eT,M,N> &mat, const BinaryExpr<eeT,OpAdd,ExprT1,ExprT2> &expr)
        : Op(mat, expr)
    {
    }

};

template <typename ExprT1, typename ExprT2>
const BinaryExpr<typename prio_type<typename element_type<ExprT1>::type,
                        typename element_type<ExprT2>::type>::result_type,
                        OpAdd, ExprT1, ExprT2> operator+(const MatrixBase<ExprT1> &expr1, const MatrixBase<ExprT2> &expr2)
{
    return OpAdd::add(expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm

#endif

