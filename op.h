#ifndef NARUTOACM_OPFUNCS_H_
#define NARUTOACM_OPFUNCS_H_

#include <iostream>
#include <algorithm>
#include <assert.h>
#include "matshell.h"
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
    template <typename eT, int M, int N>
    friend class MatrixShell;

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
                std::cout << "add a block with size " << blocksize_ << std::endl;
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
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const EntityMatrixBase<eT2,ExprT2> &expr)
    {
        return (void *)(mat.Data()) == (void *)(expr.Data());
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr)
    {
        return is_ref_to_mat(mat, expr.expr1_) || is_ref_to_mat(mat, expr.expr2_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static bool is_ref_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,OpSMul,ExprT2,ExprT3> &expr)
    {
        return is_ref_to_mat(mat, expr.expr1_);
    }

    // 0: no ref, 1: ref but ordered to mat, 2: ref and not ordered to mat
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const EntityMatrixBase<eT2,ExprT2> &expr)
    {
        return (is_ref_to_mat(mat, expr) ? (sizeof(eT1)<=sizeof(eT2)?1:2) : 0);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr)
    {
        int ar = relation_to_mat(mat, expr.expr1_);
        return (ar == 2 ? ar : std::max(ar, relation_to_mat(mat, expr.expr2_)));
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,OpSMul,ExprT2,ExprT3> &expr)
    {
        return relation_to_mat(mat, expr.expr1_);
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename ExprT2, typename ExprT3>
    static int relation_to_mat(const EntityMatrixBase<eT1,ExprT1> &mat, const BinaryExpr<eT2,M,N,OpMul,ExprT2,ExprT3> &expr)
    {
        return (is_ref_to_mat(mat, expr) ? 2 : 0);
    }
    

    // relation < 0 means that we don't know the relation between res and expr
    template <typename eT, typename ExprT1, typename ExprT2>
    static void eval(EntityMatrixBase<eT,ExprT1> &res, const EntityMatrixBase<eT,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || rows<ExprT2>::value == 0 || rows<ExprT1>::value == rows<ExprT2>::value) 
            && (cols<ExprT1>::value == 0 || cols<ExprT2>::value == 0 || cols<ExprT1>::value == cols<ExprT2>::value)>();
        
        relation = (relation >= 0 ? relation : relation_to_mat(res, expr));
        if (relation == 0)
        {
            res.SetShape(expr.Rows(), expr.Cols());
            memcpy(res.Data(), expr.Data(), expr.Elems() * sizeof(eT));
        }
    }
    template <typename eT1, typename eT2, typename ExprT1, typename ExprT2>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const EntityMatrixBase<eT2,ExprT2> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || rows<ExprT2>::value == 0 || rows<ExprT1>::value == rows<ExprT2>::value) 
            && (cols<ExprT1>::value == 0 || cols<ExprT2>::value == 0 || cols<ExprT1>::value == cols<ExprT2>::value)>();

        relation = (relation >= 0 ? relation : relation_to_mat(res, expr));

        if (relation == 2)
        {
            auto block = cache.Acquire();
            MatrixShell<eT2,rows<ExprT2>::value,cols<ExprT2>::value> ms((eT2 *)(block->ptr), block->size/sizeof(eT2));
            eval(ms, expr, cache, 0);
            eval(res, ms, cache, 0);
        }
        else
        {
            res.SetShape(expr.Rows(), expr.Cols());
            for (int i = 0; i < expr.Elems(); ++i)
                res(i) = expr(i);
        }
    }
    template <typename eT1, typename ExprT1, typename eT2, int M, int N, typename opT, typename ExprT2, typename ExprT3>
    static void eval(EntityMatrixBase<eT1,ExprT1> &res, const BinaryExpr<eT2,M,N,opT,ExprT2,ExprT3> &expr, Cache &cache, int relation = -1)
    {
        assert<(rows<ExprT1>::value == 0 || M == 0 || rows<ExprT1>::value == M) 
        && (cols<ExprT1>::value == 0 || N == 0 || cols<ExprT1>::value == N)>();
        
        opT::eval(res, expr, cache, relation);
    }

    template <typename eT, typename ExprT>
    static size_t get_cache_block_size(const EntityMatrixBase<eT,ExprT> &expr)
    {
        return sizeof(eT) * expr.Elems();
    }
    template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
    static size_t get_cache_block_size(const BinaryExpr<eT,M,N,opT,ExprT1,ExprT2> &expr)
    {
        return std::max(sizeof(eT) * expr.Elems(), std::max(get_cache_block_size(expr.expr1_), get_cache_block_size(expr.expr2_)));
    }

protected:
    template <typename eT1, typename ExprT1, typename eT2, typename ExprT2>
    Op(EntityMatrixBase<eT1,ExprT1> &mat, const MatrixBase<eT2,ExprT2> &expr)
        : cache(get_cache_block_size(expr.Derived()))
    {
        eval(mat, expr.Derived(), cache, -1);
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

