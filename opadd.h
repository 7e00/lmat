#ifndef NARUTOACM_OPADD_H_
#define NARUTOACM_OPADD_H_

#include <assert.h>
#include <string.h>
#include "matshell.h"
#include "matrix.h"
#include "binaryexpr.h"
#include "op.h"

namespace narutoacm
{

using meta::types::prio_type;
using meta::element_type;
using meta::rows;
using meta::cols;
using meta::helper::max;
using meta::helper::is_same_type;
using meta::operator_num;

class OpAdd : public Op
{
    friend class Op;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
    friend class BinaryExprHelper;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;

protected:

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int rows(const BinaryExpr<eT,M,N,OpAdd,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int cols(const BinaryExpr<eT,M,N,OpAdd,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Cols();
    }

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,OpAdd,ExprT1,ExprT2> &expr, int idx)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(idx)) + T(expr.expr2_(idx));
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,OpAdd,ExprT1,ExprT2> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(r, c)) + T(expr.expr2_(r, c));
    }

    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void add(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));

        if (relation1 == 0 && relation2 == 0)
        {
            if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
            {
                Op::eval(res, expr1.Derived(), cache, relation1);
                for (int i = 0; i < expr2.Elems(); ++i)
                    res(i) += expr2(i);
            }
            else if (is_same_type<eT,typename prio_type<eT,eT2>::result_type>::result)
            {
                Op::eval(res, expr2.Derived(), cache, relation2);
                for (int i = 0; i < expr1.Elems(); ++i)
                    res(i) += expr1(i);
            }
            else
            {
                typedef typename prio_type<eT1,eT2>::result_type T;
                res.SetShape(expr1.Rows(), expr1.Cols());
                for (int i = 0; i < expr1.Elems(); ++i)
                {
                    res(i) = T(expr1(i)) + T(expr2(i));
                }
            }
        }
        else if (relation1 == 2 || relation2 == 2)
        {
            typedef typename prio_type<eT1,eT2>::result_type T;
            auto block = cache.Acquire();
            MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            add(ms, expr1.Derived(), expr2.Derived(), cache, 0, 0);
            Op::eval(res, ms, cache, 0);
        }
        else
        {
            typedef typename prio_type<typename prio_type<eT,eT1>::result_type,eT2>::result_type T;
            res.SetShape(expr1.Rows(), expr1.Cols());
            for (int i = 0; i < expr1.Elems(); ++i)
                res(i) = T(expr1(i)) + T(expr2(i));
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void add(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const SimpleMatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));

        if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result
            && relation1 == 0 && relation2 == 0)
        {
            Op::eval(res, expr1.Derived(), cache, relation1);
            add(res, res, expr2.Derived(), cache, 1, relation2);
        }
        else
        {
            typedef typename prio_type<eT1,eT2>::result_type T;
            auto block = cache.Acquire();
            MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            add(ms, expr1.Derived(), expr2.Derived(), cache, 0, 0);
            Op::eval(res, ms, cache, 0);
        }
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void add(EntityMatrixBase<eT,ExprT> &res, const SimpleMatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        add(res, expr2.Derived(), expr1.Derived(), cache, relation2, relation1);
    }
    template <typename eT, typename ExprT, typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static void add(EntityMatrixBase<eT,ExprT> &res, const MatrixBase<eT1,ExprT1> &expr1, const MatrixBase<eT2,ExprT2> &expr2, Cache &cache, int relation1 = -1, int relation2 = -1)
    {
        relation1 = (relation1 >= 0 ? relation1 : relation_to_mat(res, expr1.Derived()));
        relation2 = (relation2 >= 0 ? relation2 : relation_to_mat(res, expr2.Derived()));
        
        if (relation1 == 0 && relation2 == 0)
        {
            if (operator_num<OpAdd,ExprT1>::value >= operator_num<OpAdd,ExprT2>::value)
            {
                Op::eval(res, expr1.Derived(), cache, 0);
                auto block = cache.Acquire();
                MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
                Op::eval(ms, expr2.Derived(), cache, 0);
                add(res, res, ms, cache, 1, 0);
            }
            else
            {
                Op::eval(res, expr2.Derived(), cache, 0);
                auto block = cache.Acquire();
                MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms((eT1 *)(block->ptr), block->size/sizeof(eT1));
                Op::eval(ms, expr1.Derived(), cache, 0);
                add(res, ms, res, cache, 0, 1);
            }
        }
        else if (relation2 == 0 && is_same_type<eT,typename prio_type<eT,eT2>::result_type>::result)
        {
            auto block = cache.Acquire();
            MatrixShell<eT1,meta::rows<ExprT1>::value,meta::cols<ExprT1>::value> ms((eT1 *)(block->ptr), block->size/sizeof(eT1));
            Op::eval(ms, expr1.Derived(), cache, 0);
            Op::eval(res, expr2.Derived(), cache, 0);
            add(res, ms, res, cache, 0, 1);
        }
        else if (relation1 == 0 && is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
        {
            auto block = cache.Acquire();
            MatrixShell<eT2,meta::rows<ExprT2>::value,meta::cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
            Op::eval(ms, expr2.Derived(), cache, 0);
            Op::eval(res, expr1.Derived(), cache, 0);
            add(res, res, ms, cache, 1, 0);
        }
        else
        {
            typedef typename prio_type<eT1,eT2>::result_type T;
            auto block = cache.Acquire();
            MatrixShell<T,meta::rows<ExprT>::value,meta::cols<ExprT>::value> ms((T *)(block->ptr), block->size/sizeof(T));
            add(ms, expr1, expr2.Derived(), cache, 0, 0);
            Op::eval(res, ms, cache, 0);
        }
    }

    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,OpAdd,ExprT2,ExprT3> &expr, Cache &cache, int relation = -1)
    {
        relation = (relation == 0 ? 0 : -1);
        add(res, expr.expr1_, expr.expr2_, cache, relation, relation);
    }
    
protected:
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    OpAdd(EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,OpAdd,ExprT2,ExprT3> &expr)
        : Op(mat, expr)
    {
    }

};

template <typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
inline const BinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result,
                        OpAdd, ExprT1, ExprT2> operator+(const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    assert(expr1.Rows() == expr2.Rows() && expr1.Cols() == expr2.Cols());
    return GenBinaryExpr<typename prio_type<expreT1,expreT2>::result_type,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result,
                        OpAdd>(expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm

#endif

