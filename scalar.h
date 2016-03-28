#ifndef NARUTOACM_SCALAR_H_
#define NARUTOACM_SCALAR_H_

namespace narutoacm
{

template<typename eT>
class Scalar
{
public:
    inline Scalar(eT val = 0)
        :val_(val)
    {
    }
    inline Scalar(const Scalar &sc)
        :val_(sc.val_)
    {
    }
    template<typename eeT>
    inline explicit Scalar(const Scalar<eeT> &sc)
        :val_((eeT)sc)
    {
    }
    inline Scalar &operator=(const Scalar &sc)
    {
        if (this != &sc)
            val_ = sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator=(const Scalar<eeT> &sc)
    {
        val_ = (eeT)sc;
        return *this;
    }

    inline Scalar &operator+=(const Scalar &sc)
    {
        val_ += sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator+=(const Scalar<eeT> &sc)
    {
        val_ += (eeT)sc;
        return *this;
    }
    inline Scalar &operator-=(const Scalar &sc)
    {
        val_ -= sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator-=(const Scalar<eeT> &sc)
    {
        val_ -= (eeT)sc;
        return *this;
    }
    inline Scalar &operator*=(const Scalar &sc)
    {
        val_ *= sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator*=(const Scalar<eeT> &sc)
    {
        val_ *= (eeT)sc;
        return *this;
    }
    inline Scalar &operator/=(const Scalar &sc)
    {
        val_ /= sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator/=(const Scalar &sc)
    {
        val_ /= (eeT)sc;
        return *this;
    }
    inline Scalar &operator%=(const Scalar &sc)
    {
        val_ %= sc.val_;
        return *this;
    }
    template<typename eeT>
    inline Scalar &operator%=(const Scalar<eeT> &sc)
    {
        val_ %= (eeT)sc;
        return *this;
    }

    inline Scalar &operator++()
    {
        ++val_;
        return *this;
    }
    inline Scalar operator++(int)
    {
        Scalar tmp = *this;
        ++val_;
        return tmp;
    }
    inline Scalar &operator--()
    {
        --val_;
        return *this;
    }
    inline Scalar operator--(int)
    {
        Scalar tmp = *this;
        --val_;
        return tmp;
    }

    inline Scalar &operator&=(const Scalar &sc)
    {
        val_ &= sc.val_;
        return *this;
    }
    inline Scalar &operator|=(const Scalar &sc)
    {
        val_ |= sc.val_;
        return *this;
    }
    inline Scalar &operator^=(const Scalar &sc)
    {
        val_ ^= sc.val_;
        return *this;
    }

    inline Scalar &operator<<=(int n)
    {
        val_ <<= n;
        return *this;
    }
    inline Scalar &operator>>=(int n)
    {
        val_ >>= n;
        return *this;
    }

    inline operator const eT() const
    {
        return val_;
    }
    inline operator eT()
    {
        return val_;
    }

private:
    eT val_;
};

} // end of namespace narutoacm

#endif
