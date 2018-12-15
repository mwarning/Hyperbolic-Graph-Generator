
/*
 * Hyperbolic Graph Generator
 *
 * Rodrigo Aldecoa, Northeastern University
 * raldecoa@neu.edu
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

use log::*;
use std::process;
use std::f64;

use crate::hg_formats::*;


fn rho(alpha: f64, R: f64, r: f64) -> f64 {
  return alpha * (alpha * (r - R)).exp();
}

fn hg_heaviside(zeta: f64, R : f64, r1: f64, r2: f64, theta: f64) -> f64 {
  let x = (1.0 / zeta)
    * ((zeta * r1).cosh() * (zeta * r2).cosh()
    - (zeta * r1).sinh() * (zeta * r2).sinh() * theta.cos()).acosh();
  if x >= R {
    return 0.0;
  } else {
    return 1.0;
  }
}

fn hg_fermi_dirac_std(beta: f64, zeta: f64, R: f64, r1: f64, r2: f64, theta: f64) -> f64 {
  let x = (1.0 / zeta)
    * ((zeta * r1).cosh() * (zeta * r2).cosh()
      -(zeta * r1).sinh() * (zeta * r2).sinh() * theta.cos()).acosh();
  
  return 1.0 / (1.0 + (beta * zeta / 2.0 * (x - R)).exp());
}

fn hg_fermi_dirac_scm(beta: f64, R: f64, r1: f64, r2: f64) -> f64 {  
  return 1.0 / (1.0 + ((beta / 2.0) * (r1 + r2 - R)).exp());
}

fn hg_integral_heaviside(x: &[f64], _dim: usize, fp: &hg_f_params) -> f64 {
  return 1.0 / HG_PI
    * rho(fp.alpha, fp.R, x[0])
    * rho(fp.alpha, fp.R, x[1])
    * hg_heaviside(fp.zeta, fp.R, x[0], x[1], x[2]);
}

fn hg_integral_standard(x: &[f64], _dim: usize, fp: &hg_f_params) -> f64 {
  return 1.0 / HG_PI 
    * rho(fp.alpha, fp.R, x[0])
    * rho(fp.alpha, fp.R, x[1])
    * hg_fermi_dirac_std(fp.beta, fp.zeta, fp.R, x[0], x[1], x[2]);
}

fn hg_integral_scm(x: &[f64], _dim: usize, fp: &hg_f_params) -> f64 {
  return rho(fp.alpha, fp.R, x[0])
    * rho(fp.alpha, fp.R, x[1])
    * hg_fermi_dirac_scm(fp.eta, fp.R, x[0], x[1]);
}

pub fn hg_get_R(pg: &hg_parameters_t, p: &hg_algorithm_parameters_t) -> f64 {
  let calls = 100000; // number of integral iterations
  let mut xl = [0.0f64; 3];
  let mut xu = [0.0f64; 3];

  // hyperbolic_rgg and hyperbolic_standard integrals are 3D
  // soft_configuration_model only 2 dimensions
  let (dim, mut params)  = match pg.gtype {
    hg_graph_type::SOFT_CONFIGURATION_MODEL => {
      (2, hg_f_params::new(0.0, p.alpha, -1.0, p.eta, -1.0))
    },
    hg_graph_type::HYPERBOLIC_STANDARD => {
      xu[2] = HG_PI;
      (3, hg_f_params::new(0.0, p.alpha, pg.zeta_eta, -1.0, 1.0 / pg.temperature))
    },
    _ => {
      xu[2] = HG_PI;
      (3, hg_f_params::new(0.0, p.alpha, pg.zeta_eta, -1.0, -1.0))
    }
  };

  let cb = match pg.gtype {
    hg_graph_type::SOFT_CONFIGURATION_MODEL => hg_integral_scm,
    hg_graph_type::HYPERBOLIC_STANDARD => hg_integral_standard,
    _ => hg_integral_heaviside
  };

  rgsl::RngType::env_setup();
  let t : rgsl::RngType = rgsl::rng::default();
  let mut r = rgsl::Rng::new(&t).unwrap();
  let mut s = rgsl::MiserMonteCarlo::new(dim).unwrap();

  let n = pg.expected_n as f64;
  let k_bar = pg.expected_degree as f64;
  
  let eps = 0.01f64; // maximum error for the avg degree
  let mut low = 0f64;
  let e = if params.beta < 1.0 { 2.5 } else { 2.0 };
  let mut high = n.ln().powf(e).max(50.0);
  let mut res = 0.0;
  let mut mid = 0.0;

  for it in 0..5000 {
    // set midpoint
    mid = (high + low) / 2.0;
    xu[0] = mid;
    xu[1] = mid;
    params.R = mid; // R = mid

    // integrate
    let (_res, _err) = s.integrate(dim, |k| { cb(k, dim, &params) }, &mut xl, &mut xu, calls, &mut r).unwrap();
    res = _res;

    if res.is_nan() {
      mid *= 1.00001; 
    } else {
      // bisection
      if (n * res) < k_bar {
       high = mid;
      } else {
       low = mid;
      }
    }

    debug!("{} - {} - {}", it, n * res, params.R);
    if ((n * res - k_bar).abs() <= eps && !res.is_nan()) || (high <= f64::MIN_POSITIVE) {
      break;
    }
  }

  if res.is_nan() || ((n * res - k_bar).abs() > eps) || (high < f64::MIN_POSITIVE) {
    error!("Network cannot be generated. Try different parameters.");
    process::exit(1);
  }

  return mid;
}

// Given that |z|>1, we need some transformations
fn hypergeometric_f(_a: f64, mut b: f64, _c: f64, z: f64) -> f64 {
  if b == 1.0 {
    return -1.0 * (1.0 - z).ln() / z;
  } else {
    let w = 1.0 / (1.0 - z);
    
    // If b is an integer, the function diverges
    // We just add a small epsilon to make it work
    if ((b as i64) as f64) == b {
      let eps = 0.000001f64;
      b += eps;
    }

    let f = w * b / (b - 1.0) * rgsl::hypergeometric::hyperg_2F1(1.0, 1.0, 2.0 - b, w)
      + b * HG_PI * (1.0 - w).powf(-b) * w.powf(b) * (b * HG_PI).sin().powf(-1.0);

    return f;
  }
}

pub fn hg_get_lambda(params: &hg_parameters_t, _p: &hg_algorithm_parameters_t) -> f64 {
  let beta = 1.0 / params.temperature;
  let n = params.expected_n as f64;
  let k_bar = params.expected_degree as f64;

  let eps = 0.001; // maximum error for the avg degree
  let mut low = 1.0;
  let mut high = rgsl::DBL_MAX;
  let mut mid;

  loop {
    // set midpoint
    mid = (high + low) / 2.0;

    let res = hypergeometric_f(1.0, 1.0 / beta, 1.0 + 1.0 / beta, -mid);

    if (n * res) < k_bar {
      high = mid;
    } else {
      low = mid;
    }

    //cout << n*res << " - " << mid << endl;
    if (n * res - k_bar).abs() <= eps && !res.is_nan() {
      break;
    }
  }
  return mid;
}
