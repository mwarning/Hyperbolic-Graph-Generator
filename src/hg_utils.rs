/*
 * Hyperbolic Graph Generator
 *
 * Chiara Orsini, CAIDA, UC San Diego
 * chiara@caida.org
 *
 * Copyright (C) 2014 The Regents of the University of California.
 *
 * This file is part of the Hyperbolic Graph Generator.
 *
 * The Hyperbolic Graph Generator is free software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or  (at your option) any later version.
 *
 * The Hyperbolic Graph Generator is distributed in the hope that it will
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with the Hyperbolic Graph Generator.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

use crate::hg_formats::*;


pub fn hg_infer_hg_type(expected_gamma: f64, temperature: f64) -> HgGraphType {
  if expected_gamma < HG_INF_GAMMA { // finite gamma
    if temperature == 0.0 {
      HgGraphType::HyperbolicRgg
    } else if temperature < HG_INF_TEMPERATURE {
      HgGraphType::HyperbolicStandard
    } else {
      HgGraphType::SoftConfigurationModel
    }
  } else { // gamma = infinite
    if temperature == 0.0 {
      HgGraphType::AngularRgg
    } else if temperature < HG_INF_TEMPERATURE {
      HgGraphType::SoftRgg
    } else {
      HgGraphType::ErdosRenyi
    }
  }
}

/* How to extract from any probability distribution:
 * The probability integral transform states that if X is
 * a continuous random variable with cumulative distribution
 * function F_X, then the random variable Y=F_X(X) has a
 * uniform distribution on [0, 1].
 * The inverse probability integral transform is just the
 * inverse of this: specifically, if Y has a uniform distribution
 * on [0, 1] and if X has a cumulative distribution F_X, then
 * the cumulative distribution function of the random variable
 * F_X^{-1}(Y) is F_X .
 * REF:  http://en.wikipedia.org/wiki/Inverse_transform_sampling
 */

/* extract a radial coordinate value uniformly from an 
 * exponential distribution (0,radius)
 */
pub fn hg_uniform_radial_coordinate(radius: f64, rnd: &mut FnMut() -> f64) -> f64 {
  // pag.4, equation (7)
  // distribution: rho(r) = [ sinh(r) / (cosh(radius)-1) ]
  if radius == 0.0 {
    eprintln!("HGG Warning: Radius = 0.");
    return 0.0;
  }

  // uniformly extracted variable (y)
  let y = rnd();
  /* CDF of rho(r) is:
   * [(cosh(r)-1) / (cosh(radius)-1)]
   * then: 
   * r = arcosh[ 1 + y * [cosh(radius)-1] ] */
  return (1.0 + y * (radius.cosh() - 1.0)).acosh();
}

/* extract a radial coordinate value quasi-uniformly from an 
 * exponential distribution (0,radius). alpha defines how
 * uniform is the extraction.
 */
pub fn hg_quasi_uniform_radial_coordinate(radius: f64, alpha: f64, rnd: &mut FnMut() -> f64) -> f64 {
  // pag.6, equation (17)
  // distribution:  rho(r) = alpha * { [alpha * sinh(r)] / cosh(alpha *radius-1) }
  if (radius == 0.0) || (alpha == 0.0) {
    eprintln!("HGG Warning: Radius = 0 or alpha = 0: discontinuity.");
    return 0.0; // alpha = 0 -> limit does not exist (discontinuity)
  }
  // uniformly extracted variable (y)
  let y = rnd();
  /* CDF of rho(r) is:
   * [(cosh(alpha *r)-1) / (alpha * cosh(alpha*radius)-1)]
   * then: 
   * r = 1/alpha * arcosh[ 1 + alpha * y * [cosh(alpha*radius)-1] ] */
  return 1.0 / alpha * (1.0 +  y * ((alpha * radius).cosh() - 1.0)).acosh();
}

/* extract an angular coordinate value uniformly at random
 * from the uniform distribution [0, 2 PI] 
 */
pub fn hg_uniform_angular_coordinate(rnd: &mut FnMut() -> f64) -> f64 {
  // pag. 4 - subsection A (in text)
  // [0,1] -> [0,2pi]
  return rnd() * 2.0 * HG_PI;
}
