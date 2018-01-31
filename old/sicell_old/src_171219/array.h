// Trying to make some array class, modeled after
// the one in the channelflow source code. The channelflow array really in turn seems to be
// modeled after the Numerical Recipes vector.

#include <assert.h>
#include <iostream>

using std::cout;



template <class T> // The array will take arguments of template type T.
class Array {      // Or rather the array will contain data of type T as well as take such arguments.
public:
        //Constructors. 
        Array(int N); //Only make space for an 1d array of size N.
        Array(int N, const T& a); //Also set all the entries to constant value a.

        //Destructors
        ~Array();

        //Operators
        inline T& operator[](int i);             // ("inline" for speed-up apparently)
        inline const T& operator[](int i) const; //defines a []-operator to get the element at index i
        Array<T> & operator=(const Array<T> &rhs); //assignment

        //Utility methods
        inline int size() const; //returns the size of the array.
        void print() const; //prints the contents to terminal

        //Assignment
        void set(const int i,const T& a);
        T* dataPtr();
        
private:
        T* data_; //The array data is stored in this method.
        int N_; //The array size.
  
}; //end class Array - remember semicolon!

//------------ Constructors and destructor ------------
template<class T>
Array<T>::Array(int N) : N_(N), data_(N>0 ? new T[N] : NULL) {}

template<class T>
Array<T>::Array(int N,const T& a) : N_(N), data_(N>0 ? new T[N] : NULL)
{
        for(int i=0;i<N;i++) data_[i] = a;
}

template <class T>
Array<T>::~Array() {
        delete[] data_;
        data_=0;
}

//------------ Operators ------------
template <class T>
inline T& Array<T>::operator[](int i){
        assert(i>=0 && i < N_);
        return data_[i];
}

template <class T>
inline const T& Array<T>::operator[](int i) const{
        assert(i>=0 && i < N_);
        return data_[i];
}

template<class T>
Array<T> & Array<T>::operator=(const Array<T> &rhs)
{
        if (this != &rhs)
        {
                if(N_ != rhs.N_)
                {
                        if (data_ != NULL)
                                delete [] (data_);
                        N_ = rhs.N_;
                        data_=N_>0 ? new T[N_] : NULL;
                }
                for (int i=0;i<N_;i++)
                        data_[i] = rhs[i];
        }
        return *this;
}

//------------ Utilities ------------
template<class T>
inline int Array<T>::size() const {
        return N_;
}

template<class T>
inline void Array<T>::print() const {
        for (int i=0;i<N_;i++)
                cout << data_[i] << '\n';
}

template<class T>
void Array<T>::set(const int i, const T& a){
        data_[i] = a;
}

template<class T>
T* Array<T>::dataPtr()
{
        return &data_[0];
}
