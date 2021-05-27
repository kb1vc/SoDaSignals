#pragma once

#include <vector>
#include <memory>



namespace SoDa {
  
  template<typename T>
  class CircularBuffer {
  public:

    CircularBuffer(int buf_len, T prime_v) : 
      buf_len(buf_len), buf_ptr(new std::vector<T>(buf_len))
    {
      if(buf_len > 0) prime(prime_v); 
    }
				    
    ~CircularBuffer() {
    }

    void resize(int _buf_len, T prime_v) {
      buf_len = _buf_len;
      buf_ptr->resize(buf_len);
      end_idx = 0; 
      if(buf_len > 0) prime(prime_v);
    }
    
    void prime(T v) {
      end_idx = 0;
      for(auto & b : *buf_ptr) {
	b = v; 
      }
    }

    const T & operator[](int idx) {
      int a_idx = bumpIndex(end_idx, idx);

      return (*buf_ptr)[a_idx]; 
    }

    void push(const T & v) {
      (*buf_ptr)[end_idx] = v;
      
      end_idx = bumpIndex(end_idx, 1);
    }

    int size() {
      return buf_len; 
    }			
				

  protected:
    
    int bumpIndex(int v, int a) {
      v = v + a; 
      while (v >= buf_len) {
	v = v - buf_len; 
      }
      return v;
    }
    
    std::shared_ptr<std::vector<T>> buf_ptr;
    int end_idx;
    int buf_len;
    bool full; 
  };

  
  
}
