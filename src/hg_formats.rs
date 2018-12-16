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


// constants definition
pub const HG_PI : f64 = std::f64::consts::PI;
pub const HG_INF_TEMPERATURE : f64 = 10.0;
pub const HG_INF_GAMMA : f64 = 10.0;
pub const HG_INF_RADIUS : f64 = 1000.0;

pub type Graph = (Vec<HgCoordinateType>, Vec<HgConnectionType>);

/* Polar coordinates of a node in a
 * hyperbolic space */
pub struct HgCoordinateType {
  pub r: f64,     // radial coordinate, distance from the origin
  pub theta: f64, // angular coordinate, angular distance from a reference point
}

/* Bidirectional connection between two nodes */
pub struct HgConnectionType {
  pub id: usize,
  pub other_id: usize,
}

pub enum HgGraphType {
  HyperbolicRgg,
  HyperbolicStandard,
  SoftConfigurationModel,
  AngularRgg,
  SoftRgg,
  ErdosRenyi
}

/* Parameters describing a graph
 * generated in a hyperbolic space */
pub struct HgParametersType {
  pub gtype: HgGraphType,
  pub expected_n: usize,
  pub temperature: f64,
  pub expected_gamma: f64,
  pub expected_degree: f64,
  pub zeta_eta: f64,
  pub starting_id: usize
}

impl HgParametersType {
  pub fn new(n: usize,
             k_bar: f64,
             exp_gamma: f64,
             t: f64,
             zeta_eta: f64,
             gtype: HgGraphType) -> Self {
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
pub struct HgAlgorithmParametersType {
  pub radius: f64,
  pub alpha: f64,
  pub eta: f64,
  pub c: f64,
}

impl HgAlgorithmParametersType {
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
pub struct HgFParams {
  pub rr: f64,
  pub alpha: f64, 
  pub zeta: f64,
  pub eta: f64,
  pub beta: f64,
}

impl HgFParams {
  pub fn new(rr: f64, alpha: f64, zeta: f64, eta: f64, beta: f64) -> Self {
    Self {
      rr: rr,
      alpha: alpha,
      zeta: zeta,
      eta: eta,
      beta: beta,
    }
  }
}
