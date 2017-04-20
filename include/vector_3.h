/*s***************************************************************************
 *
 *   Copyright (c), Ilya Valuev 2005        All Rights Reserved.
 *
 *   Author	: Ilya Valuev, MIPT, Moscow, Russia
 *
 *   Project	: GridMD, ivutils
 *
 *   $Revision: 1.29 $
 *   $Date: 2013/08/16 11:32:45 $
 *   @(#) $Header: /home/plasmacvs/source_tree/ivutils/include/vector_3.h,v 1.29 2013/08/16 11:32:45 valuev Exp $
 *
 *****************************************************************************/
/*
$Source: /home/plasmacvs/source_tree/ivutils/include/vector_3.h,v $
$Revision: 1.29 $
$Author: valuev $
$Date: 2013/08/16 11:32:45 $
*/
/*s****************************************************************************
 * $Log: vector_3.h,v $
 * Revision 1.29  2013/08/16 11:32:45  valuev
 * compiled tcpengine
 *
 * Revision 1.28  2012/10/04 16:37:45  valuev
 * sync with KINTECH svn
 *
 * Revision 1.27  2012/08/22 12:34:03  valuev
 * Made compatible with icc10
 *
 * Revision 1.26  2012/08/22 08:36:45  valuev
 * *** empty log message ***
 *
 * Revision 1.25  2012/06/29 10:50:12  valuev
 * added linear constraints
 *
 * Revision 1.57  2012/05/30 22:22:58  lesha
 * *** empty log message ***
 *
 * Revision 1.56  2012/04/24 03:12:43  lesha
 * step_fun is added
 *
 * Revision 1.55  2012/04/13 12:41:24  valuev
 * moved stdlib.h back (NULL is defined there under llinux)
 *
 * Revision 1.54  2012/04/13 08:13:21  valuev
 * reverted to M_PI
 *
 * Revision 1.53  2012/04/12 20:11:24  lesha
 * documentation
 *
 * Revision 1.52  2012/02/29 23:41:06  lesha
 * compilation with vs10
 *
*******************************************************************************/

#ifndef VECTOR_3_H
#define VECTOR_3_H

/// \en @file vector_3.h \brief N-dimensional vectors and some useful functions for work with them.
/// \ru @file vector_3.h \brief N-мерные вектора и работа с ними.

#include <stdlib.h> // нужна для NULL,  RAND_MAX в randdir
#include <cmath>
#include <limits>
#include <algorithm> // min


/*
# ifndef fmod
//r деление по модулю чисел с плавающей точкой
# define fmod(a,b)  ((a)-((long)((a)/(b))*(b)))
# endif
*/
using namespace std;

#ifndef SINGLE_PRECISION
typedef double vec_type;
#else
typedef float vec_type;
#endif

///\en "infinitely large" number
///\ru "бесконечно большое" число
#define VEC_INFTY numeric_limits<vec_type>::max()

///\en "infinitely large" integer number
#define INT_INFTY numeric_limits<int>::max()

///\en how bigger is "ininitely small" number then numeric_limits<...>::epsilon().
///    We use this parameter since numeric_limits<...>::epsilon() is too small.
///\ru во сколько раз "бесконечно малое" число больше numeric_limits<...>::epsilon().
///    Этот множитель необходим потому, что numeric_limits<...>::epsilon() слишком мало.
#define MULT_EPSILON 1024

///\en "infinitely small" number
///\ru "бесконечно малое" число
#define VEC_ZERO MULT_EPSILON*numeric_limits<vec_type>::epsilon()

///\en Storage template prototype for a vector.
class storage_prototype{
  ///\en Storage class should define value_type, double used as example here.
  typedef double value_type;
  ///\en Storage class should have default constructor
  storage_prototype();

  ///\en Storage class should have copy constructor
  storage_prototype(storage_prototype &other);

  ///\en Storage class should have resize function for setting new size.
  ///    It is called for "result" vectors created by default constructor to communicate actual dimension,
  ///    implementation may be skipped if dimension is preset. New elements must be initialized by 0.
  void resize(size_t new_dim);

  ///\en Same as \ref resize() but with a value to be set for all new elements.
  void resize(size_t new_dim, const value_type& value);

  ///\en Storage class should provide function for querying vector dimension.
  size_t dim() const;
  
  ///\en Storage class should provide [] operator for randomly accessing its elements for read and write (size_t index).
  value_type& operator[](size_t i) const ;

  ///\en Storage class should provide [] operator for randomly accessing its elements for read and write (integer index).
  ///    Integer index is used for compatibility with old Vector_3 code.
  value_type& operator[](int i) const ;

  
};


///\en Array storage class for a vector followinfg \ref storage_prototype.
template <class T, size_t N>
struct array_stor_t{
  typedef T value_type;
  
  ///\en The data is contiguous array in memory.
  T data[N];

  ///\en The dimension is preset to N, this function does nothing.
  void resize(size_t new_dim){}

  void resize(size_t new_dim, const T & value){}

  inline size_t dim() const { return N; }
  
  inline T& operator[](int i) const {return (T&)data[i];}

  inline T& operator[](size_t i) const {return (T&)data[i];}

  ///\en gets address in memory for pointer association
  T *get_ptr() { return data; }

  const T *get_ptr() const { return data; }
};



///\en General N-dimensional vector of type T and some useful operations.
///    Storage class for a vector should follow \ref storage_prototype.
///    Note that for binary operators on vectors with variable dimension storages 
///    the result dimension is inherited from the argument with the minimal dimension.
///\ru N-мерный вектор типа T с некоторыми полезными операциями
template <class T, size_t N, class storage_t= array_stor_t<T,N> > 
class Vector_Nt{
public:
  storage_t v;
  
public:  
  
  typedef T value_type;

  size_t dim() const { return v.dim(); } 

  Vector_Nt(const storage_t &stor): v(stor){}

  ///\en makes all components equal to a
  ///\ru присваивает всем компонентам значение a
  Vector_Nt(const T &a=0){
    for(size_t i=0;i<v.dim();i++) v[i]=a;
  }

  inline Vector_Nt& operator=(const T &a){
    for(size_t i=0;i<v.dim();i++) v[i]=a;
    return *this;
  }

  explicit Vector_Nt(const T &a1, const T &a2){
    if(v.dim()>0)v[0]=a1;
    if(v.dim()>1)v[1]=a2;
    for(size_t i=2;i<v.dim();i++)v[i]=0;
  }

  explicit Vector_Nt(const T &a1, const T &a2, const T& a3){
    if(v.dim()>0)v[0]=a1;
    if(v.dim()>1)v[1]=a2;
    if(v.dim()>2)v[2]=a3;
    for(size_t i=3;i<v.dim();i++)v[i]=0;
  }

  ///\en construct from input iterator (or array)
  template <class A>
  Vector_Nt(const A *beg){
    for(size_t i=0;i<v.dim();i++,++beg)
      v[i]=*beg;
  }

  ///\en Construct from another vector, if dimension is not preset, inherits \a vect dimension.
  template <class T2, size_t M, class storage2_t>
  Vector_Nt(const Vector_Nt<T2,M,storage2_t> &vect){
    v.resize(vect.dim()); // has effect only if dimension is variable
    size_t dim_=min(v.dim(),vect.dim());
    for(size_t i=0;i<dim_;i++) v[i]=vect[i];
  }

  ///\en Copy from another vector.
  ///    If dimensions are different, copies  min(dim1,dim2) number of elements only.
  template <class T2, size_t M, class storage2_t>
  inline Vector_Nt& operator=(const Vector_Nt<T2,M,storage2_t> &vect){
    size_t dim_=min(v.dim(),vect.dim());
    for(size_t i=0;i<dim_;i++) v[i]=vect[i];
    return *this;
  }

  ///\en Copies vector to iterator
  ///\ru Копирует содержимое вектора в итератор
  template <class A>
  void copy_to(A *beg) const{  
    for(size_t i=0;i<v.dim();i++,++beg)
      *beg=v[i];
  }

  ///\en obtains element value
  ///\ru получение элемента 
  inline T& operator[](int i) const {return (T&)v[i];}

  inline T& operator[](size_t i) const {return (T&)v[i];}

  ///\en comparison. If the difference is less then VEC_ZERO (MULT_EPSILON*numeric_limits<T>::epsilon()) then components are assumed to be equal
  ///\ru сравнение. При отличии меньше чем на VEC_ZERO (MULT_EPSILON*numeric_limits<T>::epsilon()) компоненты считаются одинаковыми
  template <size_t M, class storage2_t>
  inline bool operator==(const Vector_Nt<T, M, storage2_t>& vect) const{
    size_t dim_=min(v.dim(),vect.dim());
    for(size_t i=0;i<dim_;i++)
      if(fabs(v[i]-vect[i])>MULT_EPSILON*numeric_limits<T>::epsilon())return false;
    return true;
  }

  template <size_t M, class storage2_t>
  inline bool operator!=(const Vector_Nt<T, M, storage2_t>& vect) const{
    return (!(this->operator==(vect)));
  }

  template <size_t M, class storage2_t>
  inline Vector_Nt operator+(const Vector_Nt<T, M, storage2_t>& vect) const{
    Vector_Nt result;
    size_t dim_=min(v.dim(),vect.dim());
    result.v.resize(dim_); // has effect only for variable dimension storage
    for (size_t i=0; i<dim_ ;i++)
      result.v[i]=v[i]+vect.v[i];
    return result;
  }

  template <size_t M, class storage2_t>
  inline Vector_Nt operator-(const Vector_Nt<T, M, storage2_t>& vect) const {
    Vector_Nt result;
    size_t dim_=min(v.dim(),vect.dim());
    result.v.resize(dim_); // has effect only for variable dimension storage
    for (size_t i=0; i<dim_ ;i++)
      result.v[i]=v[i]-vect.v[i];
    return result;
  }
 
  ///\en Scalar product
  ///\ru Скалярное произведение векторов
  template <size_t M, class storage2_t>
  inline T operator*(const Vector_Nt<T, M, storage2_t>& vect) const{
    T result=0;
    size_t dim_=min(v.dim(),vect.dim());
    for (size_t i=0; i<dim_; i++)
      result+=v[i]*vect.v[i];
    return result;
  }

  ///\en Multiplies on coefficient
  ///\ru Покомпонентное умножение на коэффициент
  inline Vector_Nt operator*(const T &coeff) const{
    Vector_Nt result;
    result.v.resize(v.dim()); // has effect only for variable dimension storage
    for (size_t i=0; i<v.dim(); i++)
      result[i]=coeff*v[i];
    return result;
  }

  ///\en Multiplies on vector by components
  ///\ru Покомпонентное умножение на вектор
  template <class Vec>
  inline Vector_Nt scale(const Vec &s) const{
    Vector_Nt result;
    size_t dim_=min(v.dim(),s.dim());
    result.v.resize(dim_); // has effect only for variable dimension storage
    for(size_t i=0;i<v.dim();i++)
      result[i]=s[i]*v[i];
    return result;
  }


  ///\en Divides on coefficient
  ///\ru Покомпонентное деление на коэффициент
  inline Vector_Nt operator/(const T &coeff) const {
    Vector_Nt result;
    result.v.resize(v.dim()); // has effect only for variable dimension storage
    for (size_t i=0; i<v.dim(); i++)
      result[i]=v[i]/coeff;
    return result;
  }

  ///\en Vector product (dim()=3 only)
  ///\ru Векторное произведение
  template <size_t M, class storage2_t>
  inline Vector_Nt operator%(const Vector_Nt<T, M, storage2_t>& vect) const{ //reserved for N specializations
    if(v.dim()==3)
      return Vector_Nt(v[1]*vect.v[2]-v[2]*vect.v[1],v[2]*vect.v[0]-v[0]*vect.v[2],v[0]*vect.v[1]-v[1]*vect.v[0]);
    return *this;
  }

  ///\en Multiplies on -1
  ///\ru Умножение вектора на -1
  inline Vector_Nt operator-() const {
    Vector_Nt result;
    result.v.resize(v.dim()); // has effect only for variable dimension storage
    for(size_t i=0;i<v.dim();i++)
      result.v[i]=-v[i];
    return result;
  }

  ///\ru Сложение с присваиванием
  template <size_t M, class storage2_t>
  inline Vector_Nt& operator+=(const Vector_Nt<T, M, storage2_t> &vect){
    size_t dim_=min(v.dim(),vect.dim());
    for (size_t i=0; i<dim_; i++)
      v[i]+=vect.v[i];
    return *this;
  }

  ///\ru Вычитание с присваиванием
  template <size_t M, class storage2_t>
  inline Vector_Nt& operator-=(const Vector_Nt<T, M, storage2_t> &vect){
    size_t dim_=min(v.dim(),vect.dim());
    for (size_t i=0; i<dim_; i++)
      v[i]-=vect.v[i];
    return *this;
  }

  ///\ru Умножение на коэффициент с присваиванием
  inline Vector_Nt& operator*=(const T &coeff){
    for (size_t i=0; i<v.dim(); i++)
      v[i]*=coeff;
    return *this;
  }

  ///\ru Деление на скаляр с присваиванием
  inline Vector_Nt& operator/=(const T &coeff){
    for (size_t i=0; i<v.dim(); i++)
      v[i]/=coeff;
    return *this;
  }

  ///\en Norm squared
  ///\ru Квадрат нормы вектора
  T norm2() const {
    T result=0;
    for (size_t i=0; i<v.dim(); i++)
      result+=v[i]*v[i];
    return result; 
  }

  ///\en Norm
  ///\ru Норма вектора
  T norm() const {
    return sqrt(norm2());
  }

  ///\en Returns norm and normalizes vector on newnorm
  ///\ru Возвращает норму и нормализует вектор на newnorm
  T normalize(T newnorm=1.){
    T norm=this->norm();
    if(norm>=MULT_EPSILON*numeric_limits<T>::epsilon()){
      T c=newnorm/norm;
      for (size_t i=0; i<v.dim(); i++)
        v[i]*=c;
    }
    return norm;
  }

  ///\en Normalizes a vector that may have infinite components
  T inormalize(T newnorm=1., T infty=VEC_INFTY){
    T result=0;
    for(size_t i=0;i<v.dim();i++){
      if(fabs(v[i])>=infty){
        if(result>=0){
          for(size_t j=0;j<i;j++)
            v[j]=0.;
          result=0.;
        }
        result-=1;
        v[i]=v[i]>0 ? 1 : -1;
      }
      else if(result>=0) //still summing the normal components
        result+=v[i]*v[i];
      else
        v[i]=0.;
    }
    if(fabs(result)<MULT_EPSILON*numeric_limits<T>::epsilon())
      return 0.;
    newnorm/=sqrt(fabs(result));
    for (size_t i=0; i<v.dim(); i++)
      v[i]*=newnorm;
    return result<0 ? infty : result;
  }

  ///\en return a vector projection on a given axis
  Vector_Nt prj(size_t i) const {
    Vector_Nt res;
    res.v.resize(v.dim(),0.); // has effect only for variable dimension storage
    res[i]=v[i];
    return res;
  }

  ///\en return a vector projection on a given vector (k*v)*k
  Vector_Nt prj(const Vector_Nt &k) const {
    return k*(k*(*this));
  }

  ///\en return a vector of length l parallel to this
  Vector_Nt unit(T l=1) const {
    Vector_Nt res(*this);
    res.normalize(l);
    return res;
  }

  ///\en nearest image distance within rectangular cell (FOR DISTANCE MEASUREMENTS)
  ///    assumes that each coordinate absolute value is in the range [0,cell[i]) 
  ///    returned vector is in the range [-cell[i]/2,cell[i]/2)
  ///    flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  ///    Note that \a cell dimension is not checked.
  ///\ru Ближайший образ в прямоугольной ячейке
  ///    Считаем, что все пространство разделено на ячейки - параллелепипеды с ребрами, параллельными
  ///    осям координат и диагональю, заданной вектором rcell, причем начало координат является 
  ///    центром одной из ячеек. Если *this находится центральной ячейке, возвращается копия *this.\n
  ///    Иначе, если *this находится в кубе 3*3 ячейки с центром в начале координат, то создает образ 
  ///    *this в центральной ячейке.\n
  ///    Иначе, возвращает неопределенное значение.  
  template <size_t M, class storage2_t>
  Vector_Nt rcell1(const Vector_Nt<T, M, storage2_t> &cell,int flags=0xffff) const{
    Vector_Nt ret(*this);
    for(size_t i=0;i<v.dim();i++){
      if(flags&(0x1<<i)){
        if(v[i]>cell[i]/2){
          ret[i]-=cell[i];
        }
        else if(v[i]<-cell[i]/2){
          ret[i]+=cell[i];
        }
      }
    }
    return ret;
  }

  ///\en reduction to elementary cell [0, cell[i]) (FOR REDUCTION TO ELEMENTARY CELL)
  ///    flags indicate the periodicity in specific directions: 0x1 for X, 0x2 for Y, 0x4 for Z
  ///    Note that \a cell dimension is not checked.
  ///\ru Почти то же, что и rcell1, но без ограничения на положение *this и с другой системой ячеек.
  ///    В начале координат находится не центр ячейки, а ее угол. Может работать медленнее из-за наличия 
  ///    операции деления по модулю с плавающей точкой
  template <size_t M, class storage2_t>
  Vector_Nt rcell(const Vector_Nt<T, M, storage2_t> &cell, int flags=0xffff) const {
    Vector_Nt ret(*this);
    for (size_t i=0, flag=1; i<v.dim(); i++, flag<<=1) {
      if(flags&flag){
        ret.v[i]=fmod(v[i],cell.v[i]);
        if(ret[i]<0)ret.v[i]+=cell.v[i];
      }
    }
    return ret;
  }

  ///\en the same as rcell, but start point of zero cell is p1
  ///    Note that \a p1,  \a cell dimensions are not checked.
  ///\ru то же самое, что и rcell, только начало нулевой ячейки имеет координаты p1
  template <size_t M, class storage2_t, size_t P, class storage3_t>
  Vector_Nt rpcell(const Vector_Nt<T, M, storage2_t> &p1, const Vector_Nt<T, P, storage3_t> &cell, int flags=0xfff) const {
    Vector_Nt ret(*this);
    for (size_t i=0, flag=1; i<v.dim(); i++, flag<<=1) {
      if(flags&flag){
        if (ret[i]<p1[i] || ret[i]>p1[i]+cell[i]) {
          ret[i]=fmod(v[i]-p1[i],cell[i])+p1[i];
          if (ret[i]<p1[i]) ret[i]+=cell[i];
        }
      }
    }
    return ret;
  }
  
  ///\en returns maximal vector component and its index
  ///\ru Возвращает максимальную компоненту вектора и ее индекс в ind 
  T maxcoord(int *ind=NULL) const {
    if(!v.dim())
      return 0.;
    int im=0;
    T vv=v[0];
    for (size_t i=1; i<v.dim(); i++) {
      if(v[i]>vv){
        im=(int)i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }

  ///\en returns minimal vector component and its index
  ///\ru Возвращает минимальную компоненту вектора и ее индекс в ind 
  T mincoord(int *ind=NULL) const {
    if(!v.dim())
      return 0.;
    int im=0;
    T vv=v[0];
    for (int i=1; i<v.dim(); i++) {
      if(v[i]<vv){
        im=(int)i;
        vv=v[i];
      }
    }
    if(ind)*ind=im;
    return vv;
  }

  ///\en returns the coord having maximal absolute value
  T maxabscoord(int *ind=NULL) const {
    if(!v.dim())
      return 0.;
    int im=0;
    T vv=fabs(v[0]);
    for (size_t i=1; i<v.dim(); i++) {
      if(fabs(v[i])>vv){
        im=i;
        vv=fabs(v[i]);
      }
    }
    if(ind)*ind=im;
    return v[im];
  }

  ///\en returns the coord having minimal absolute value
  T minabscoord(int *ind=NULL) const {
    if(!v.dim())
      return 0.;
    int im=0;
    T vv=fabs(v[0]);
    for (size_t i=1; i<v.dim(); i++) {
      if(fabs(v[i])<vv){
        im=(int)i;
        vv=fabs(v[i]);
      }
    }
    if(ind)*ind=im;
    return v[im];
  }

  ///\en returns true if the vector has infinite components
  bool infinite(T infty=VEC_INFTY) const {
    for(size_t i=0;i<v.dim();i++){
      if(fabs(v[i])>=infty)
        return true;
    }
    return false;
  }

  T *get_ptr() { return v.get_ptr(); }

  const T *get_ptr() const { return v.get_ptr(); }


/*
  //r Выводит вектор в поток вывода по умолчанию, в формате (x,y,z)\\n
  void print() const{
    cout<< "(";
    for(size_t i=0;i<v.dim();i++){
      cout<< v[i];
      if(i!=v.dim()-1)
        cout<< ", ";
    }
    cout<< ")\n";
  }*/
};


template<class T, size_t N, class storage_t>
Vector_Nt<T, N, storage_t> operator*(const T &coeff,const Vector_Nt<T, N, storage_t> &vec){
  return vec*coeff;
}

//template<class T, class T2, size_t N, class storage_t>
//Vector_Nt<T, N, storage_t> operator*(const T2 &coeff,const Vector_Nt<T, N, storage_t> &vec){
//  return vec*coeff;
//}


// old Vector_3 compatibility typedefs and functions
typedef Vector_Nt<int,2> iVector_2;
typedef Vector_Nt<int,3> iVector_3;
typedef Vector_Nt<vec_type, 2> Vector_2;
typedef Vector_Nt<vec_type, 3> Vector_3;
typedef Vector_3 *Vector_3P; 
typedef Vector_2 *Vector_2P;


template <size_t N> 
class  Vector_N: public Vector_Nt<vec_type, N>::type{
};

///\en returns polygon area based on vectors vect1 and vect2
///\ru возвращает площадь параллелограмма, построенного на векторах vect1 и vect2
inline vec_type vec_area(const Vector_2 &vect1, const Vector_2 &vect2) {
  return fabs(vect1[0]*vect2[1]-vect1[1]*vect2[0]);
};

inline vec_type vec_area(const Vector_3 &vect1, const Vector_3 &vect2) {
  return (vect1%vect2).norm();
};



///\en finds the maximum distance between vector pairs
///\ru Находит максимальное расстояние между векторами va1[i], va2[i], i=1..n
///    \param va1 - массив Vector_3[n]
///    \param n - длина массивов va1 и va2
vec_type dist_max(Vector_3 *va1,Vector_3 *va2,int n);

///\en finds average distance between vector pairs
///\ru Находит среднее расстояние между векторами va1[i], va2[i], i=1..n
vec_type dist_av(Vector_3 *va1,Vector_3 *va2,int n);

///\en finds the average difference norm between two vector sets of the same length
///    optionally gives the indices for maximal and minimal difference
///    va2 can be NULL, then the norm of va1 is used
///\ru Находит среднее расстояние между va1[i] и va2[i], а также, по желанию, индексы, на которых достигается min и max расстояние
vec_type diff_av(Vector_3 *va1,Vector_3 *va2,int n, int *minind=0, int *maxind=0);

///\en finds suitable perpendicular to a vector
///\ru Находит перпендикуляр к вектору vAB
Vector_3 FindPerp(const Vector_3 &vAB);

///\en Returns the average (center) vector of the vector array
///    and cooordinates of a minimal box in which it is contained
Vector_3 GetScope(const Vector_3 *varr,long n,Vector_3* box_min,Vector_3* box_max);

///\en the same with long index array
Vector_3 GetIScope(const Vector_3 *varr,long *indarr,long n,Vector_3* box_min,Vector_3* box_max);

///\en the same with int index array
Vector_3 GetIScopei(const Vector_3 *varr,int *indarr,int n,Vector_3* box_min,Vector_3* box_max);

// neue Funktionen

///\en clears vector array with optional integer index
///\ru Очистка массива векторов, с поддержкой индексирования 
///    В данном Vector_3 vec[] обнуляет n координат. Если ind==NULL, то 
///    очищает первые n элементов. Если ind!=NULL, то для i=0..n-1
///    очищает vec[ind[i]]
///    См. \ref indexed_calculations.
void clear_vecarri(int n,Vector_3 *vec, int *ind=0);

///\en reflects the vector ini+dir*t+0.5*force*t^2 to be inside a box limited by 0 and box sizes
///    changes dir according to the final state
///    fills crossed dir with bit flags corresponding directions along which the walls were crossed
Vector_3 Reflect(Vector_3& ini, double t,Vector_3 &dir, double *box, int flag=0x7, const Vector_3 &force=Vector_3()); 

///\en returns random unit vector uniformely distributed in space (?? check this)
Vector_3 randdir();



///\en Calculates extent of the vector container.
///    \return the center of the vector set, optionally
///    (if arguments are not NULL) fills the bounding box in \a box_min, \a box_max.
template<class vec_inp_it>
Vector_3 get_extent(vec_inp_it beg,vec_inp_it end, Vector_3* box_min=NULL,Vector_3* box_max=NULL){
  if(beg==end)
    return Vector_3();
  Vector_3 center(*beg++);
  Vector_3 cube1(center), cube2(center);
  size_t n=1;
  for(;beg!=end;++beg){
    Vector_3 vec=*beg;
    center+=vec;
    for(size_t j=0;j<3;j++){
      if(cube1[j]>vec[j])
        cube1[j]=vec[j];
     if(cube2[j]<vec[j])
        cube2[j]=vec[j];
    }
    n++;
  }
  if(box_min)
    *box_min=cube1;
  if(box_max)
    *box_max=cube2;
  return center/n;
}


///\en Performs a step of the Stabilized Gramm-Schidt orthonormalization algorithm. 
///    Given a set of orthogonal unit vectors defined by [orth_beg, orth_end)
///    orthonormalizes inp_vec with respect to this set (removes all projections to the set).
///    The result is recorded to out_vec, which is normalized.
///    \return the norm of out_vec before normalization. If return value is <=\ref VEC_ZERO,
///    then inp_vec may be considered as linearly dependent of the supplied set.
template<class inp_it, class vector_tt>
typename vector_tt::value_type gramm_schmidt_project(inp_it orth_set_beg, inp_it orth_set_end, const vector_tt &inp_vec, vector_tt &out_vec, typename vector_tt::value_type new_norm = 1.){
  out_vec = inp_vec;
  for(;orth_set_beg!=orth_set_end; ++orth_set_beg){
    out_vec -= ( (vector_tt)(*orth_set_beg) *out_vec)*( (vector_tt)(*orth_set_beg) ); 
  }
  if(new_norm>0.)
    return out_vec.normalize(new_norm);
  else
    return 0.;
}


# endif // VECTOR_3_H

