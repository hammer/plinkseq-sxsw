#ifndef __PSEQ_MATRIX_H__
#define __PSEQ_MATRIX_H__

#include <vector>
#include "helper.h"

namespace Data { 

  template<class T> class Matrix;

  template<class T = double> class Vector {

    public:

    Vector() { } 
    Vector(int i) { resize(i); }
    Vector(const std::vector<T> & x ) { data = x; }
    
    void clear() { data.clear(); }
    void resize(const int i) { data.resize(i); }
    void resize(const int i, const T & t) { data.resize(i,t); }
    void push_back( const T & t ) { data.push_back(t); }
    int size() const { return data.size(); }
    int dim1() const { return data.size(); }
    
    T  operator[] ( const unsigned int i) const { return data[i]; } 
    T & operator[] ( const unsigned int i) { return data[i]; } 
    
    T operator() (const unsigned int i ) const { return data[i]; }
    T & operator() (const unsigned int i ) { return data[i]; }
    
    Vector<T> operator*( const Matrix<T> & rhs ) const;
    Vector<T> operator+( const Vector<T> & rhs ) const;
    Vector<T> operator-( const Vector<T> & rhs ) const;
    
    private:
    
    std::vector<T> data;

  };

  

  template<class T = double> class Matrix {
    
    public:

    // row access

    struct Row
    { 
      Row( Matrix & m , int i ) : mat(m), row(i) { }
      T & operator[](const int j) { return mat(row,j); }     
      //T operator[](const int j) const { return mat(row,j); }     
      private:
      int row;
      Matrix & mat;      
    };
    
    struct ConstRow
    { 
      ConstRow( const Matrix & m , int i ) : mat(m) , row(i) { } 
      T operator[](const int j) const { return mat(row,j); }           
      private:
      int row;
      const Matrix & mat;      
    };
    
    Matrix() { clear(); } 
    Matrix(const int r, const int c) { resize(r,c); }
    
    T operator() (const unsigned int i, const unsigned int j ) const { return data[j][i]; }
    T & operator() (const unsigned int i, const unsigned int j ) { return data[j][i]; }
    
    Row operator[] ( const unsigned int i) { return Row(*this,i); }
    ConstRow operator[] ( const unsigned int i) const { return ConstRow(*this,i); }
    
    Vector<T> row( const int r ) 
    { 
      Vector<T> d( ncol );
      for (int c=0; c<ncol; c++) d[c] = (*this)(r,c);
      return d;
    } 

    Vector<T> col( const int c ) const { return data[c]; } 
    Vector<T> & col( const int c ) { return data[c]; } 

    void add_col( const Vector<T> & r ) { data.push_back(r); ++ncol; }
    void add_col( const std::vector<T> & r ) { data.push_back( Vector<T>(r) ); ++ncol; }
    
    void cbind( const Data::Matrix<T> & rhs )
    {
      if ( nrow != rhs.dim1() ) 
	Helper::halt( "cbind() for matrices with unequal number of rows" );
      for (int c=0; c<rhs.dim2(); c++)
	add_col( rhs.col(c) );
    }
    
    void add_row( const Vector<T> & r ) 
    { 
      if ( r.size() != ncol ) 
	{ 
	  if ( nrow == 0 ) resize(0,r.size());
	  else { plog.warn("bad row addition"); return; }
	}

      for( int i=0; i<ncol; i++ ) data[i].push_back( r[i] );
      ++nrow;
    }

    void add_row( const std::vector<T> & r ) 
    { 
      if ( r.size() != ncol ) 
	{
	  if ( nrow == 0 ) resize(0,r.size());
	  else { plog.warn("bad row addition"); return; }
	}

      for( int i=0; i<ncol; i++ ) data[i].push_back( r[i] );
      ++nrow;
    }
    
    void set_row_mask( int r , const bool b = true ) 
    { 
      if ( r >= 0 && r < nrow ) row_mask[r] = b;
    }
    
    bool masked( const int r )
    {
      if ( ncol == 0 ) return false; // no data
      if ( r >= 0 && r < nrow ) return row_mask[r];
      return true; // mask out-of-range items
    }

    void clear() { data.clear(); row_mask.clear(); nrow = ncol = 0; }
    
    void resize(const int r, const int c) 
    { 
      nrow = r;
      ncol = c;
      row_mask.resize( nrow , false ); // masked-out
      data.resize(c); 
      for (int j=0; j<c; j++) data[j].resize( nrow );
    }
    
    int dim1() const { return nrow; }
    int dim2() const { return ncol; }


    // op. overloading for common matrix operations
    
    Matrix<T> operator*( const Data::Matrix<T> & rhs ) const;
    Vector<T> operator*( const Data::Vector<T> & rhs ) const;
    
    private:
    
    std::vector< Vector<T> > data;
    std::vector<bool> row_mask;
    int nrow ;
    int ncol ;
  };


  
  

}


#endif
