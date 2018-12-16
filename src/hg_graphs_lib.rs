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

use crate::hg_gen_algorithms::*;
use crate::hg_formats::*;
use crate::hg_utils::*;


pub fn hg_hyperbolic_distance(
        gtype: hg_graph_type,
        node1: &hg_coordinate_t,
        node2: &hg_coordinate_t,
        zeta_eta: f64) -> f64 {
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }

  match gtype {
    hg_graph_type::HYPERBOLIC_RGG | hg_graph_type::HYPERBOLIC_STANDARD =>
      hg_hyperbolic_distance_hyperbolic_rgg_standard(zeta_eta, node1, node2),
    hg_graph_type::SOFT_CONFIGURATION_MODEL =>
      hg_hyperbolic_distance_scm(node1, node2),
    hg_graph_type::ANGULAR_RGG | hg_graph_type::SOFT_RGG =>
      hg_hyperbolic_distance_angular_soft_rgg(node1, node2),
    hg_graph_type::ERDOS_RENYI =>
      hg_hyperbolic_distance_er(node1, node2)
  }
}

pub fn hg_graph_generator(
        n: usize,
        k_bar: f64, 
        exp_gamma: f64,
        t : f64, 
        zeta: f64,
        seed: u32) -> Result<Graph, &'static str> {
  // Define a uniform random number distribution which produces
  // "double" values between 0 and 1 (0 inclusive, 1 exclusive).
  // Equivalent to the boost::mt19937 generator used in the cpp version.
  rgsl::RngType::env_setup();
  let rnd_type = rgsl::rng::algorithms::mt19937();
  let mut rnd = rgsl::Rng::new(&rnd_type).unwrap();

  // set seed
  rnd.set(seed as usize);

  let mut rnd_01 = || { rnd.uniform() };
  let gt = hg_infer_hg_type(exp_gamma, t);

  match gt {
    hg_graph_type::HYPERBOLIC_RGG =>
      hg_hyperbolic_rgg(n, &mut rnd_01, k_bar, exp_gamma, zeta),
    hg_graph_type::HYPERBOLIC_STANDARD =>
      hg_hyperbolic_standard(n, &mut rnd_01, k_bar, exp_gamma, t, zeta),
    hg_graph_type::SOFT_CONFIGURATION_MODEL =>
      hg_soft_configuration_model(n, &mut rnd_01, k_bar, exp_gamma, zeta),
    hg_graph_type::ANGULAR_RGG =>
      hg_angular_rgg(n, &mut rnd_01, k_bar, zeta),
    hg_graph_type::SOFT_RGG =>
      hg_soft_rgg(n, &mut rnd_01, k_bar, t, zeta),
    hg_graph_type::ERDOS_RENYI =>
      hg_erdos_renyi(n, &mut rnd_01, k_bar, zeta),
  }
}
