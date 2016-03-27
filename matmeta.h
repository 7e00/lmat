#ifndef NARUTOACM_MATMETA_H_
#define NARUTOACM_MATMETA_H_

namespace narutoacm
{
template <typename eT>
class Scalar;
template <typename eT, int M, int N>
class Matrix;
template <typename eT, int M, int N>
class MatrixShell;
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
class BinaryExpr;
template <typename eT, int M, int N, typename opT, typename ExprT>
class UnaryExpr;


class MatOp;
class MatOpAdd;
class MatOpSub;
class MatOpEMul; // element-wise multiplication
class MatOpEDiv;
class MatOpSMul; // matrix * scalar
class MatOpMul; // matrix multiplication
class MatOpTran; // transpose


typedef double f64;
typedef float f32;
typedef long long s64;
typedef unsigned long long u64;
typedef int s32;
typedef unsigned int u32;
typedef short s16;
typedef unsigned short u16;
typedef signed char s8;
typedef unsigned char u8;


namespace meta
{
namespace helper
{
template <bool>
struct assert;
template<>
struct assert<true> {};

template <int a, int b>
struct max
{
    enum { result = (a > b ? a : b) };
};
template <int a, int b>
struct min
{
    enum { result = (a < b ? a : b) };
};

template <typename eT1,typename eT2>
struct is_same_type
{
    enum { result = false };
};
template <typename eT>
struct is_same_type<eT,eT>
{
    enum { result = true };
};

template <bool, typename eT1, typename eT2>
struct select_type;
template <typename eT1, typename eT2>
struct select_type<true,eT1,eT2>
{
    typedef eT1 result_type;
};
template <typename eT1, typename eT2>
struct select_type<false,eT1,eT2>
{
    typedef eT2 result_type;
};

} // end of namespace narutoacm::meta::helper

namespace types
{
template <typename eT1, typename eT2>
struct prio_type { typedef eT1 result_type; };
template<typename eT> struct prio_type<eT,f64> { typedef f64 result_type; };
template<> struct prio_type<s64,f32> { typedef f32 result_type; };
template<> struct prio_type<s32,f32> { typedef f32 result_type; };
template<> struct prio_type<s16,f32> { typedef f32 result_type; };
template<> struct prio_type<s8,f32> { typedef f32 result_type; };
template<> struct prio_type<u64,f32> { typedef f32 result_type; };
template<> struct prio_type<u32,f32> { typedef f32 result_type; };
template<> struct prio_type<u16,f32> { typedef f32 result_type; };
template<> struct prio_type<u8,f32> { typedef f32 result_type; };
template<> struct prio_type<u64,s64> { typedef s64 result_type; };
template<> struct prio_type<s32,s64> { typedef s64 result_type; };
template<> struct prio_type<s16,s64> { typedef s64 result_type; };
template<> struct prio_type<s8,s64> { typedef s64 result_type; };
template<> struct prio_type<u32,s64> { typedef s64 result_type; };
template<> struct prio_type<u16,s64> { typedef s64 result_type; };
template<> struct prio_type<u8,s64> { typedef s64 result_type; };
template<> struct prio_type<s32,u64> { typedef s64 result_type; };
template<> struct prio_type<s16,u64> { typedef s64 result_type; };
template<> struct prio_type<s8,u64> { typedef s64 result_type; };
template<> struct prio_type<u64,s32> { typedef s64 result_type; };
template<> struct prio_type<u64,s16> { typedef s64 result_type; };
template<> struct prio_type<u64,s8> { typedef s64 result_type; };
template<> struct prio_type<u32,u64> { typedef u64 result_type; };
template<> struct prio_type<u16,u64> { typedef u64 result_type; };
template<> struct prio_type<u8,u64> { typedef u64 result_type; };
template<> struct prio_type<u32,s32> { typedef s32 result_type; };
template<> struct prio_type<u16,s32> { typedef s32 result_type; };
template<> struct prio_type<u8,s32> { typedef s32 result_type; };
template<> struct prio_type<s16,s32> { typedef s32 result_type; };
template<> struct prio_type<s8,s32> { typedef s32 result_type; };
template<> struct prio_type<u32,s16> { typedef s32 result_type; };
template<> struct prio_type<u32,s8> { typedef s32 result_type; };
template<> struct prio_type<s16,u32> { typedef s32 result_type; };
template<> struct prio_type<s8,u32> { typedef s32 result_type; };
template<> struct prio_type<u16,u32> { typedef u32 result_type; };
template<> struct prio_type<u8,u32> { typedef u32 result_type; };
template<> struct prio_type<u16,s16> { typedef s16 result_type; };
template<> struct prio_type<u8,s16> { typedef s16 result_type; };
template<> struct prio_type<s8,s16> { typedef s16 result_type; };
template<> struct prio_type<u16,s8> { typedef s16 result_type; };
template<> struct prio_type<s8,u16> { typedef s16 result_type; };
template<> struct prio_type<u8,u16> { typedef u16 result_type; };
template<> struct prio_type<u8,s8> { typedef s8 result_type; };

using helper::is_same_type;
using helper::select_type;

template <typename eT1, typename eT2>
struct inferior_type
{
    typedef typename select_type<
        is_same_type<eT1,typename prio_type<eT1,eT2>::result_type>::result,
        eT2,eT1>::result_type result_type;
};
} // end of namespace narutoacm::meta::types

using helper::assert;
using helper::max;
using helper::min;

template <typename ExprT>
struct element_type;
template <typename eT>
struct element_type<Scalar<eT>>
{
    typedef eT type;
};
template <typename eT, int M, int N>
struct element_type<Matrix<eT,M,N>>
{
    typedef eT type;
};
template <typename eT, int M, int N>
struct element_type<MatrixShell<eT,M,N>>
{
    typedef eT type;
};
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
struct element_type<BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
    typedef eT type;
};
template <typename eT, int M, int N, typename opT, typename ExprT>
struct element_type<UnaryExpr<eT,M,N,opT,ExprT>>
{
    typedef eT type;
};

template <typename ExprT>
struct rows;
template <typename ExprT>
struct cols;
template <typename eT, int M, int N>
struct rows<Matrix<eT,M,N>>
{
    static const int value = M;
};
template <typename eT, int M, int N>
struct cols<Matrix<eT,M,N>>
{
    static const int value = N;
};
template <typename eT, int M, int N>
struct rows<MatrixShell<eT,M,N>>
{
    static const int value = M;
};
template <typename eT, int M, int N>
struct cols<MatrixShell<eT,M,N>>
{
    static const int value = N;
};
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
struct rows<BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
    static const int value = M;
};
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
struct cols<BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
    static const int value = N;
};
template <typename eT, int M, int N, typename opT, typename ExprT>
struct rows<UnaryExpr<eT,M,N,opT,ExprT>>
{
    static const int value = M;
};
template <typename eT, int M, int N, typename opT, typename ExprT>
struct cols<UnaryExpr<eT,M,N,opT,ExprT>>
{
    static const int value = N;
};

template <typename ExprT>
struct matrix_expr_type;
template <typename eT, int M, int N>
struct matrix_expr_type<Matrix<eT,M,N>>
{
    enum { value = 2 };
};
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
struct matrix_expr_type<BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>> : public 
                                                    assert<(rows<ExprT1>::value == 0
                                                    || rows<ExprT2>::value == 0
                                                    || rows<ExprT1>::value == rows<ExprT2>::value)
                                                    && (M == 0
                                                    || max<rows<ExprT1>::value,rows<ExprT2>::value>::result == 0
                                                    || M == max<rows<ExprT1>::value,rows<ExprT2>::value>::result)
                                                    && (cols<ExprT1>::value == 0
                                                    || cols<ExprT2>::value == 0
                                                    || cols<ExprT1>::value == cols<ExprT2>::value)
                                                    && (N == 0
                                                    || max<cols<ExprT1>::value,cols<ExprT2>::value>::result == 0
                                                    || N == max<cols<ExprT1>::value,cols<ExprT2>::value>::result)>
{
    enum { value = min<min<1,matrix_expr_type<ExprT1>::value>::result, matrix_expr_type<ExprT2>::value>::result };
};
template <typename eT, int M, int N,  typename ExprT1, typename SeT1>
struct matrix_expr_type<BinaryExpr<eT,M,N,MatOpSMul,ExprT1,Scalar<SeT1>>> : public 
                                                    assert<(M == 0
                                                    || rows<ExprT1>::value == 0
                                                    || M == rows<ExprT1>::value)
                                                    && (N == 0
                                                    || cols<ExprT1>::value == 0
                                                    || N == cols<ExprT1>::value)>
{
    enum { value = min<1, matrix_expr_type<ExprT1>::value>::result };
};
template <typename eT, int M, int N, typename ExprT1, typename ExprT2>
struct matrix_expr_type<BinaryExpr<eT,M,N,MatOpMul,ExprT1,ExprT2>> : public
                                                    assert<(cols<ExprT1>::value == 0
                                                    || rows<ExprT2>::value == 0
                                                    || cols<ExprT1>::value == rows<ExprT2>::value)
                                                    && (M == 0
                                                    || rows<ExprT1>::value == 0
                                                    || M == rows<ExprT1>::value)
                                                    && (N == 0
                                                    || cols<ExprT2>::value == 0
                                                    || N == cols<ExprT2>::value)>
{
    enum { value = 0 };
};
template <typename eT, int M, int N, typename opT, typename ExprT>
struct matrix_expr_type<UnaryExpr<eT,M,N,opT,ExprT>> : public
                                                    assert<(M == 0
                                                    || rows<ExprT>::value == 0
                                                    || M == rows<ExprT>::value)
                                                    && (N == 0
                                                    || cols<ExprT>::value == 0
                                                    || N == cols<ExprT>::value)>
{
    enum { value = 0 };
};
template <typename eT, int M, int N, typename ExprT>
struct matrix_expr_type<UnaryExpr<eT,M,N,MatOpTran,ExprT>> : public
                                                    assert<(M == 0
                                                    || cols<ExprT>::value == 0
                                                    || M == cols<ExprT>::value)
                                                    && (N == 0
                                                    || rows<ExprT>::value == 0
                                                    || N == rows<ExprT>::value)>
{
    enum { value = min<1, matrix_expr_type<ExprT>::value>::result };
};

template <typename opT, typename ExprT>
struct operator_num;
template <typename opT, typename eT, int M, int N>
struct operator_num<opT,Matrix<eT,M,N>>
{
    static const int value = 1;
};
template <typename opT, typename eT, int M, int N>
struct operator_num<opT,MatrixShell<eT,M,N>>
{
    static const int value = 1;
};
template <typename opT, typename eT>
struct operator_num<opT,Scalar<eT>>
{
    static const int value = 1;
};
template <typename opT, typename eT, int M, int N, typename ExprT1, typename ExprT2>
struct operator_num<opT,BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
    static const int value = operator_num<opT,ExprT1>::value + operator_num<opT,ExprT2>::value;
};
template <typename eT, int M, int N, typename opT, typename ExprT1, typename ExprT2>
struct operator_num<MatOp,BinaryExpr<eT,M,N,opT,ExprT1,ExprT2>>
{
    static const int value = operator_num<MatOp,ExprT1>::value + operator_num<MatOp,ExprT2>::value;
};
template <typename opT, typename eT, int M, int N, typename ExprT>
struct operator_num<opT,UnaryExpr<eT,M,N,opT,ExprT>>
{
    static const int value = operator_num<opT,ExprT>::value;
};
template <typename eT, int M, int N, typename opT, typename ExprT>
struct operator_num<MatOp,UnaryExpr<eT,M,N,opT,ExprT>>
{
    static const int value = operator_num<MatOp,ExprT>::value;
};

} // end of namespace narutoacm::meta

} // end of namespace narutoacm

#endif

