#ifndef NARUTOACM_OPEDIV_H_
#define NARUTOACM_OPEDIV_H_

#include <string.h>
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

class OpEDiv : public Op
{
    friend class Op;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2, int type>
    friend class BinaryExprHelper;
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    friend class BinaryExpr;

protected:

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int rows(const BinaryExpr<eT,M,N,OpEDiv,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Rows();
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static int cols(const BinaryExpr<eT,M,N,OpEDiv,ExprT1,ExprT2> &expr)
    {
        return expr.expr1_.Cols();
    }

    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,OpEDiv,ExprT1,ExprT2> &expr, int idx)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(idx)) / T(expr.expr2_(idx));
    }
    template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
    static eT eval(const BinaryExpr<eT,M,N,OpEDiv,ExprT1,ExprT2> &expr, int r, int c)
    {
        typedef typename prio_type<typename prio_type<typename element_type<ExprT1>::type, eT>::result_type, typename element_type<ExprT2>::type>::result_type T;
        return T(expr.expr1_(r, c)) / T(expr.expr2_(r, c));
    }

    template <typename eT, typename eT1, int M1, int N1, typename eT2, int M2, int N2>
    static void ediv(eT *res, int len, const Matrix<eT1,M1,N1> &expr1, const Matrix<eT2,M2,N2> &expr2, Cache &cache)
    {
        if (is_same_type<eT,typename prio_type<eT,eT1>::result_type>::result)
        {
            Op::eval(res, len, expr1, cache);
            auto ptr = expr2.Data();
            for (int i = 0; i < expr2.Elems(); ++i)
                res[i] /= ptr[i];
        }
        else if (is_same_type<eT,typename prio_type<eT,eT2>::result_type>::result)
        {
            Op::eval(res, len, expr2, cache);
            auto ptr = expr1.Data();
            for (int i = 0; i < expr1.Elems(); ++i)
                res[i] = ptr[i] / res[i];
        }
        else
        {
            typedef typename prio_type<eT1,eT2>::result_type T;
            for (int i = 0; i < expr1.Elems(); ++i)
                res[i] = T(expr1(i)) / T(expr2(i));
        }
    }

    template <typename eT, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(eT *res, int len, const SimpleMatrixBase<expreT1,ExprT1> &expr1, const SimpleMatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        if (is_same_type<eT,typename prio_type<eT,expreT1>::result_type>::result)
        {
            Op::eval(res, len, mat1, cache);
            for (int i = 0; i < mat2.Elems(); ++i)
                res[i] /= mat2(i);
        }
        else if (is_same_type<eT,typename prio_type<eT,expreT2>::result_type>::result)
        {
            Op::eval(res, len, mat2, cache);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = mat1(i) / res[i];
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = T(mat1(i)) / T(mat2(i));
        }
    }

    template <typename eT, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(eT *res, int len, const MatrixBase<expreT1,ExprT1> &expr1, const SimpleMatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        if (is_same_type<eT,typename prio_type<eT,expreT1>::result_type>::result)
        {
            Op::eval(res, len, mat1, cache);
            for (int i = 0; i < mat2.Elems(); ++i)
                res[i] /= mat2(i);
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            auto ptr = (T *)(block->ptr);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = ptr[i];
        }
    }

    template <typename eT, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(eT *res, int len, const SimpleMatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        if (is_same_type<eT,typename prio_type<eT,expreT2>::result_type>::result)
        {
            Op::eval(res, len, mat2, cache);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = mat1(i) / res[i];
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            auto ptr = (T *)(block->ptr);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = ptr[i];
        }
    }

    template <typename eT, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(eT *res, int len, const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();
        
        if (is_same_type<eT,typename prio_type<eT,expreT1>::result_type>::result)
        {
            Op::eval(res, len, mat1, cache);
            auto block = cache.Acquire();
            Op::eval((expreT2 *)(block->ptr), block->size / sizeof(expreT2), mat2, cache);
            expreT2 *ptr = (expreT2 *)(block->ptr);
            for (int i = 0; i < mat2.Elems(); ++i)
                res[i] /= ptr[i];
        }
        else if (is_same_type<eT,typename prio_type<eT,expreT2>::result_type>::result)
        {
            Op::eval(res, len, mat2, cache);
            auto block = cache.Acquire();
            Op::eval((expreT1 *)(block->ptr), block->size / sizeof(expreT1), mat1, cache);
            auto *ptr = (expreT1 *)(block->ptr);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = ptr[i] / res[i];
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            auto ptr = (T *)(block->ptr);
            for (int i = 0; i < mat1.Elems(); ++i)
                res[i] = ptr[i];
        }
    }

    template <typename eT, int M, int N, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(Matrix<eT,M,N> &res, const SimpleMatrixBase<expreT1,ExprT1> &expr1, const SimpleMatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();

        int ar = relation_to_mat(res, mat1);
        int br = relation_to_mat(res, mat2);

        if (ar == 0 && br == 0)
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            ediv(res.Data(), res.Capacity(), mat1, mat2, cache);
        }
        else if (ar == 2 || br == 2)
        {
            typedef typename prio_type<typename prio_type<eT,expreT1>::result_type,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto ptr = (T *)(block->ptr);
            if (is_same_type<eT,T>::result)
            {
                memcpy(res.Data(), ptr, res.Elems()*sizeof(eT));
            }
            else
            {
                auto resptr = res.Data();
                for (int i = 0; i < mat1.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
        else // res's shape must already ok
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto resptr = res.Data();
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            for (int i = 0; i < mat1.Elems(); ++i)
                resptr[i] = T(mat1(i)) / T(mat2(i));
        }
    }

    template <typename eT, int M, int N, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(Matrix<eT,M,N> &res, const MatrixBase<expreT1,ExprT1> &expr1, const SimpleMatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();

        int ar = relation_to_mat(res, mat1);
        int br = relation_to_mat(res, mat2);

        if (ar == 0 && br == 0)
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            ediv(res.Data(), res.Capacity(), mat1, mat2, cache);
        }
        else if (br == 2)
        {
            typedef typename prio_type<typename prio_type<eT,expreT1>::result_type,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto ptr = (T *)(block->ptr);
            if (is_same_type<eT,T>::result)
            {
                memcpy(res.Data(), ptr, res.Elems()*sizeof(eT));
            }
            else
            {
                auto resptr = res.Data();
                for (int i = 0; i < mat1.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            auto block = cache.Acquire();
            Op::eval((expreT1 *)(block->ptr), block->size / sizeof(expreT1), mat1, cache);
            res.SetShape(mat2.Rows(), mat2.Cols());
            auto ptr = (expreT1 *)(block->ptr);
            auto resptr = res.Data();
            for (int i = 0; i < mat2.Elems(); ++i)
                resptr[i] = T(ptr[i]) / T(mat2(i));
        }
    }

    template <typename eT, int M, int N, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(Matrix<eT,M,N> &res, const SimpleMatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();

        int ar = relation_to_mat(res, mat1);
        int br = relation_to_mat(res, mat2);

        if (ar == 0 && br == 0)
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            ediv(res.Data(), res.Capacity(), mat1, mat2, cache);
        }
        else if (ar == 2)
        {
            typedef typename prio_type<typename prio_type<eT,expreT1>::result_type,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto ptr = (T *)(block->ptr);
            if (is_same_type<eT,T>::result)
            {
                memcpy(res.Data(), ptr, res.Elems()*sizeof(eT));
            }
            else
            {
                auto resptr = res.Data();
                for (int i = 0; i < mat1.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
        else
        {
            typedef typename prio_type<expreT1,expreT2>::result_type T;
            auto block = cache.Acquire();
            Op::eval((expreT2 *)(block->ptr), block->size / sizeof(expreT2), mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto ptr = (expreT2 *)(block->ptr);
            auto resptr = res.Data();
            for (int i = 0; i < mat2.Elems(); ++i)
                resptr[i] = T(mat1(i)) / T(ptr[i]);
        }
    }

    template <typename eT, int M, int N, typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
    static void ediv(Matrix<eT,M,N> &res, const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2, Cache &cache)
    {
        auto &mat1 = expr1.Derived();
        auto &mat2 = expr2.Derived();

        int ar = relation_to_mat(res, mat1);
        int br = relation_to_mat(res, mat2);

        if (ar == 0 && br == 0)
        {
            res.SetShape(mat1.Rows(), mat1.Cols());
            ediv(res.Data(), res.Capacity(), mat1, mat2, cache);
        }
        else if (ar == 0)
        {
            if (is_same_type<eT,typename prio_type<eT,expreT1>::result_type>::result)
            {
                auto block = cache.Acquire();
                Op::eval((expreT2 *)(block->ptr), block->size / sizeof(expreT2), mat2, cache);
                res.SetShape(mat1.Rows(), mat1.Cols());
                Op::eval(res.Data(), res.Capacity(), mat1, cache);
                auto ptr = (expreT2 *)(block->ptr);
                auto resptr = res.Data();
                for (int i = 0; i < res.Elems(); ++i)
                    resptr[i] /= ptr[i];
            }
            else
            {
                typedef typename prio_type<expreT1,expreT2>::result_type T;
                auto block = cache.Acquire();
                ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
                res.SetShape(mat1.Rows(), mat1.Cols());
                auto ptr = (T *)(block->ptr);
                auto resptr = res.Data();
                for (int i = 0; i < mat1.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
        else if (br == 0)
        {
            if (is_same_type<eT,typename prio_type<eT,expreT2>::result_type>::result)
            {
                auto block = cache.Acquire();
                Op::eval((expreT1 *)(block->ptr), block->size / sizeof(expreT1), mat1, cache);
                res.SetShape(mat2.Rows(), mat2.Cols());
                Op::eval(res.Data(), res.Capacity(), mat2, cache);
                auto ptr = (expreT1 *)(block->ptr);
                auto resptr = res.Data();
                for (int i = 0; i < res.Elems(); ++i)
                    resptr[i] = ptr[i] / resptr[i];
            }
            else
            {
                typedef typename prio_type<expreT1,expreT2>::result_type T;
                auto block = cache.Acquire();
                ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
                res.SetShape(mat1.Rows(), mat1.Cols());
                auto ptr = (T *)(block->ptr);
                auto resptr = res.Data();
                for (int i = 0; i < mat2.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
        else
        {
            typedef typename prio_type<typename prio_type<eT,expreT1>::result_type,expreT2>::result_type T;
            auto block = cache.Acquire();
            ediv((T *)(block->ptr), block->size / sizeof(T), mat1, mat2, cache);
            res.SetShape(mat1.Rows(), mat1.Cols());
            auto ptr = (T *)(block->ptr);
            if (is_same_type<eT,T>::result)
            {
                memcpy(res.Data(), ptr, res.Elems()*sizeof(eT));
            }
            else
            {
                auto resptr = res.Data();
                for (int i = 0; i < mat1.Elems(); ++i)
                    resptr[i] = ptr[i];
            }
        }
    }

    template <typename eT1, typename eT2, int M, int N, typename ExprT1, typename ExprT2>
    static void eval(eT1 *res, int len, const BinaryExpr<eT2,M,N,OpEDiv,ExprT1,ExprT2> &expr, Cache &cache)
    {
        ediv(res, len, expr.expr1_, expr.expr2_, cache);
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename ExprT1, typename ExprT2>
    static void eval(Matrix<eT1,M1,N1> &res, const BinaryExpr<eT2,M2,N2,OpEDiv,ExprT1,ExprT2> &expr, Cache &cache)
    {
        ediv(res, expr.expr1_, expr.expr2_, cache);
    }
    
protected:
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename ExprT1, typename ExprT2>
    OpEDiv(Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,OpEDiv,ExprT1,ExprT2> &expr)
        : Op(mat, expr)
    {
    }

};

template <typename expreT1, typename ExprT1, typename expreT2, typename ExprT2>
inline const BinaryExpr<double,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result,
                        OpEDiv, ExprT1, ExprT2> ElementDivide(const MatrixBase<expreT1,ExprT1> &expr1, const MatrixBase<expreT2,ExprT2> &expr2)
{
    assert(expr1.Rows() == expr2.Rows() && expr1.Cols() == expr2.Cols());
    return GenBinaryExpr<double,
                        max<meta::rows<ExprT1>::value,meta::rows<ExprT2>::value>::result,
                        max<meta::cols<ExprT1>::value,meta::cols<ExprT2>::value>::result,
                        OpEDiv>(expr1.Derived(), expr2.Derived());
}

} // end of namespace narutoacm

#endif

