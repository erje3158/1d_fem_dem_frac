#ifndef SHAPE_H
#define SHAPE_H

#include "realtypes.h"
#include "vec_frac.h"

namespace dem_frac {
  
  class shape{ // abstract class
  public:
    virtual vec  getCenter() const =0;
    virtual REAL getVolume() const =0;
    virtual vec  randomPoint() const =0;
    virtual void print() const =0;
    virtual ~shape() {} // virtual destructor
  };
  
} // namespace dem

#endif
