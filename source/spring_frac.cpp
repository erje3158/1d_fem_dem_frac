#include "spring_frac.h"
#include "parameter_frac.h"

namespace dem_frac{

spring::spring(particle_frac &ptcl1, particle_frac &ptcl2, REAL modulus)
  :p1(ptcl1), p2(ptcl2), Young(modulus) {
  init(p1, p2);
}

spring::spring(std::vector<particle_frac*> &ParticleVec, int id1, int id2, REAL modulus)
  :p1(*ParticleVec[id1]), p2(*ParticleVec[id2]), Young(modulus) {
  init(p1, p2);
}

vec spring::getDeformation() {
  vec dir =  p2.getCurrPosition() - p1.getCurrPosition();
  REAL length = vfabs( dir );
  return (length - length0) * normalize(dir);
}

void spring::applyForce() {// p1 adds this force, p2 subtracts this force
  vec dir =  p2.getCurrPosition() - p1.getCurrPosition();
  REAL length = vfabs( dir );
  if (length - length0 > EPS ) {
    vec force = (length - length0) * normalize(dir);
    force *= ks;
    p1.addForce(force);
    p2.addForce(-force);
  }
}

}
