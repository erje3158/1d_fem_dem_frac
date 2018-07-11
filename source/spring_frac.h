#ifndef SPRING_H
#define SPRING_H

#include "realtypes.h"
#include "particle_frac.h"
#include "vec_frac.h"
#include <cassert>

namespace dem_frac{

class spring {

 public:
  spring(particle_frac &p1, particle_frac &p2, REAL young);
  spring(std::vector<particle_frac*> &ParticleVec, int id1, int id2, REAL young);

  REAL getLength0() const {return length0;}
  REAL getLength() const {return vfabs( p2.getCurrPosition() - p1.getCurrPosition() );}
  vec  getDeformation(); 
  void applyForce();
  int getParticleId1() const {return p1.getID();}
  int getParticleId2() const {return p2.getID();}

 private:
  particle_frac &p1;
  particle_frac &p2;
  REAL Young;   // Young's modulus
  REAL ks;      // stiffness
  REAL length0; // equilibrium length

  void init(particle_frac &p1, particle_frac &p2) {
    length0 = vfabs( p2.getCurrPosition() - p1.getCurrPosition() );
    REAL radius = p1.getMaxRadius();
    assert (radius == p2.getMaxRadius() );
    ks = Young * 4 * radius * radius / length0;
  } 
};



}

#endif
