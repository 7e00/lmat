#ifndef NARUTOACM_SHAREDPOOL_H_
#define NARUTOACM_SHAREDPOOL_H_

#include <stack>
#include <memory>

namespace narutoacm
{
namespace util
{
template <class T>
class SharedPool
{
private:
    class External_Deleter
    {
    public:
        External_Deleter()
        {}
        explicit External_Deleter(std::weak_ptr<SharedPool<T>* > pool)
            : pool_(pool) {}

        void operator()(T* ptr)
        {
            if (auto pool_ptr = pool_.lock())
            {
                try {
                    (*pool_ptr.get())->Add(std::unique_ptr<T>{ptr});
                    return;
                } catch(...) {}
            }
            std::default_delete<T>{}(ptr);
        }
    private:
        std::weak_ptr<SharedPool<T>* > pool_;
    };

public:
    using ptr_type = std::unique_ptr<T, External_Deleter >;

    SharedPool() : this_ptr_(new SharedPool<T>*(this)) {}
    virtual ~SharedPool(){}

    void Add(std::unique_ptr<T> t)
    {
        pool_.push(std::move(t));
    }

    ptr_type Acquire()
    {
        assert(!pool_.empty());
        ptr_type tmp(pool_.top().release(),
                    External_Deleter{std::weak_ptr<SharedPool<T>*>{this_ptr_}});
        pool_.pop();
        return std::move(tmp);
    }

    bool Empty() const
    {
        return pool_.empty();
    }

    size_t Size() const
    {
        return pool_.size();
    }

private:
    std::shared_ptr<SharedPool<T>* > this_ptr_;
    std::stack<std::unique_ptr<T> > pool_;
};

} // end of namespace narutoacm::util
} // end of namespace narutoacm

#endif

