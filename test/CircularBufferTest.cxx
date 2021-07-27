#include "../include/CircularBuffer.hxx"
#include <iostream>

template<typename T>
void dumpCBuf(SoDa::CircularBuffer<T> b) {
  for(int i = 0; i < b.size(); i++) {
    std::cout << i << ": " << b[i] << "\n";
  }
}

int main(int argc, char * argv[])
{
  // create a circular buffer of chars.
  SoDa::CircularBuffer<char> chbuf(8, '\000');

  for(char c = 'a'; c < 'i'; c++) {
    chbuf.push(c);
  }
  
  dumpCBuf(chbuf);
  
  std::cout << "Push i\n";
  chbuf.push('i');
  dumpCBuf(chbuf);

  std::cout << "Push j to q\n";
  for(char c = 'j'; c <= 'q'; c++) {
    chbuf.push(c);
  }
  dumpCBuf(chbuf);
}


