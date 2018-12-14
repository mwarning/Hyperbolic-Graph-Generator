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

use std::fmt;
use log::*;

// constants definition
pub const HG_PI : f64 = std::f64::consts::PI; // 3.14159265359
pub const HG_INF_TEMPERATURE : f64 = 10.0;
pub const HG_INF_GAMMA : f64 = 10.0;
pub const HG_INF_RADIUS : f64 = 1000.0;

pub type Graph = (Vec<hg_coordinate_t>, Vec<hg_connection_t>);

/* Polar coordinates of a node in a
 * hyperbolic space */
pub struct hg_coordinate_t {
  pub r: f64,     // radial coordinate, distance from the origin
  pub theta: f64, // angular coordinate, angular distance from a reference point
}

/* Bidirectional connection between two nodes */
pub struct hg_connection_t {
  pub id: usize,
  pub other_id: usize,
}

pub enum hg_graph_type {
  HYPERBOLIC_RGG,
  HYPERBOLIC_STANDARD,
  SOFT_CONFIGURATION_MODEL,
  ANGULAR_RGG, 
  SOFT_RGG,
  ERDOS_RENYI
}

/* Parameters describing a graph
 * generated in a hyperbolic space */
pub struct hg_parameters_t {
  pub gtype: hg_graph_type, // type => gtype
  pub expected_n: usize,
  pub temperature: f64,
  pub expected_gamma: f64,
  pub expected_degree: f64,
  pub zeta_eta: f64,
  pub starting_id: usize
}

impl hg_parameters_t {
  pub fn new(n: usize,
             k_bar: f64,
             exp_gamma: f64,
             t: f64,
             zeta_eta: f64,
             gtype: hg_graph_type) -> Self {
    debug!("\tGraph initialization");
    /* initialize the graph structure with the 
     * parameters provided in input by the user */

    // Init random generator
    //HG_Random::init(seed); // TODO: move out and create hg_graph_t::new()

    Self {
      gtype: gtype,
      expected_n: n,
      temperature: t,
      expected_gamma: exp_gamma,
      expected_degree: k_bar,
      zeta_eta: zeta_eta,
      starting_id: 1,
    }
  }
}

/* Graph generation internal parameters */
pub struct hg_algorithm_parameters_t {
  pub radius: f64,
  pub alpha: f64,
  pub eta: f64,
  pub c: f64,
}

impl fmt::Debug for hg_algorithm_parameters_t {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
      write!(f, " hg_algorithm_parameters_t {{ radius: {}, alpha: {}, eta: {}, c: {}}}",
        self.radius, self.alpha, self.eta, self.c)
    }
}

impl hg_algorithm_parameters_t {
  pub fn new() -> Self {
    Self {
      // -1.0 => not relevant for current model
      radius: -1.0,
      alpha: -1.0,
      eta: -1.0,
      c: -1.0
    }
  }
}

/* Structures used for numerical *
 * integration. Not all fields   *
 * are used every time           */
pub struct hg_f_params {
  pub R: f64,
  pub alpha: f64, 
  pub zeta: f64,
  pub eta: f64,
  pub beta: f64,
}

impl hg_f_params {
  pub fn new(R: f64, alpha: f64, zeta: f64, eta: f64, beta: f64) -> Self {
    Self {
      R: R,
      alpha: alpha,
      zeta: zeta,
      eta: eta,
      beta: beta,
    }
  }
}
