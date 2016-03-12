#ifndef NARUTOACM_OPFUNCS_H_
#define NARUTOACM_OPFUNCS_H_

#include <iostream>
#include <algorithm>
#include <assert.h>
#include "matrix.h"
#include "binaryexpr.h"
#include "sharedpool.h"

namespace narutoacm
{
using meta::element_type;
using meta::helper::assert;
using meta::rows;
using meta::cols;
using util::SharedPool;

class Op
{
    template <typename eT, int M, int N>
    friend class Matrix;

protected:
    class Cache
    {
    public:
        struct Block
        {
            Block(unsigned long sz)
                : ptr(nullptr), size(sz)
            {
                ptr = new unsigned char[sz];
            }
            ~Block()
            {
                delete [] ptr;
            }
            unsigned char *ptr;
            unsigned long size;
        };
        Cache(unsigned long blocksize)
            : blocksize_(blocksize)
        {
        }
        
        using ptr_type = SharedPool<Block>::ptr_type;
        
        ptr_type Acquire()
        {
            if (pool_.Empty())
            {
                //std::cout << "add a block with size " << blocksize_ << std::endl;
                pool_.Add(std::unique_ptr<Block>(new Block(blocksize_)));
            }
            return pool_.Acquire();
        }
        
        size_t Size() const
        {
            return pool_.Size();
        }
        
    private:
        unsigned long blocksize_;
        SharedPool<Block> pool_;
    };

protected:
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    static bool is_ref_to_mat(const Matrix<eT,M,N> &mat, const MatrixBase<expreT,ExprT> &expr)
    {
        return false;
    }
    template <typename eT, int M, int N>
    static bool is_ref_to_mat(const Matrix<eT,M,N> &mat, const Matrix<eT,M,N> &expr)
    {
        return (&mat == &expr);
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename opT, typename ExprT1, typename ExprT2>
    static bool is_ref_to_mat(const Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,opT,ExprT1,ExprT2> &expr)
    {
        return is_ref_to_mat(mat, expr.expr1_) || is_ref_to_mat(mat, expr.expr2_);
    }

    // 0: no ref, 1: ref but ordered to mat, 2: ref and not ordered to mat
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    static int relation_to_mat(const Matrix<eT,M,N> &mat, const MatrixBase<expreT,ExprT> &expr)
    {
        return 0;
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2>
    static int relation_to_mat(const Matrix<eT1,M1,N1> &mat, const Matrix<eT2,M2,N2> &expr)
    {
        return ((void *)(&mat) == (void *)(&expr) ? 1 : 0);
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename opT, typename ExprT1, typename ExprT2>
    static int relation_to_mat(const Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,opT,ExprT1,ExprT2> &expr)
    {
        int ar = relation_to_mat(mat, expr.expr1_);
        return (ar == 2 ? ar : std::max(ar, relation_to_mat(mat, expr.expr2_)));
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename ExprT1, typename SeT2>
    static int relation_to_mat(const Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,OpSMul,ExprT1,Scalar<SeT2>> &expr)
    {
        return relation_to_mat(mat, expr.expr1_);
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename ExprT1, typename ExprT2>
    static int relation_to_mat(const Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,OpMul,ExprT1,ExprT2> &expr)
    {
        return (is_ref_to_mat(mat, expr) ? 2 : 0);
    }

    template <typename eT, int M, int N>
    static void eval(eT *data, int len, const Matrix<eT,M,N> &expr, Cache &cache)
    {
        assert(len >= expr.Elems());
        memcpy(data, expr.Data(), expr.Elems() * sizeof(eT));
    }
    template <typename eT, typename eeT, int M, int N>
    static void eval(eT *data, int len, const Matrix<eeT,M,N> &expr, Cache &cache)
    {
        assert(len >= expr.Elems());
        for (int i = 0; i < expr.Elems(); ++i)
            data[i] = expr(i);
    }
    template <typename eT, typename eeT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    static void eval(eT *res, int len, const BinaryExpr<eeT,M,N,opT,ExprT1,ExprT2> &expr, Cache &cache)
    {
        assert(len >= expr.Elems());
        opT::eval(res, len, expr.Derived(), cache);
    }
    template <typename eT, int M, int N>
    static void eval(Matrix<eT,M,N> &res, const Matrix<eT,M,N> &expr, Cache &cache)
    {
        if (&res != &expr)
        {
            res.SetShape(expr.Rows(), expr.Cols());
            memcpy(res.Data(), expr.Data(), expr.Elems() * sizeof(eT));
        }
    }
    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2>
    static void eval(Matrix<eT1,M1,N1> &res, const Matrix<eT2,M2,N2> &expr, Cache &cache)
    {
        assert<(M1 == 0 || M2 == 0 || M1 == M2) && (N1 == 0 || N2 == 0 || N1 == N2)>();
        res.SetShape(expr.Rows(), expr.Cols());
        eval(res.Data(), res.Capacity(), expr, cache);
    }

    template <typename eT1, int M1, int N1, typename eT2, int M2, int N2, typename opT, typename ExprT1, typename ExprT2>
    static void eval(Matrix<eT1,M1,N1> &mat, const BinaryExpr<eT2,M2,N2,opT,ExprT1,ExprT2> &expr, Cache &cache)
    {
        assert<(M1 == 0 || M2 == 0 || M1 == M2) && (N1 == 0 || N2 == 0 || N1 == N2)>();
        opT::eval(mat, expr, cache);
    }

    template <typename eT, int M, int N>
    static size_t get_cache_block_size(const Matrix<eT,M,N> &expr)
    {
        return sizeof(eT) * expr.Elems();
    }

    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    static size_t get_cache_block_size(const BinaryExpr<eT,M,N,opT,ExprT1,ExprT2> &expr)
    {
        return std::max(sizeof(eT) * expr.Elems(), std::max(get_cache_block_size(expr.expr1_), get_cache_block_size(expr.expr2_)));
    }

protected:
    template <typename eT, int M, int N, typename expreT, typename ExprT>
    Op(Matrix<eT,M,N> &mat, const MatrixBase<expreT,ExprT> &expr)
        : cache(get_cache_block_size(expr.Derived()))
    {
        eval(mat, expr.Derived(), cache);
    }

private:
    Cache cache;
};

template <typename eT, int M, int N>
template <typename expreT, typename ExprT>
Matrix<eT,M,N>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
{
    Op op(*this, matexpr.Derived());
}
template <typename eT, int M, int N>
template <typename expreT, typename ExprT>
Matrix<eT,M,N> &Matrix<eT,M,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT>
template <typename expreT, typename ExprT>
Matrix<eT,0,0>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), row_(0), col_(0)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT>
template <typename expreT, typename ExprT>
Matrix<eT,0,0> &Matrix<eT,0,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT, int M>
template <typename expreT, typename ExprT>
Matrix<eT,M,0>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), col_(0)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT, int M>
template <typename expreT, typename ExprT>
Matrix<eT,M,0> &Matrix<eT,M,0>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT, int N>
template <typename expreT, typename ExprT>
Matrix<eT,0,N>::Matrix(const MatrixBase<expreT,ExprT> &matexpr)
    : data_(nullptr), capacity_(0), row_(0)
{
    Op op(*this, matexpr.Derived());
}

template <typename eT, int N>
template <typename expreT, typename ExprT>
Matrix<eT,0,N> &Matrix<eT,0,N>::operator=(const MatrixBase<expreT,ExprT> &matexpr)
{
    Op op(*this, matexpr.Derived());
}

} // end of namespace narutoacm



#endif

