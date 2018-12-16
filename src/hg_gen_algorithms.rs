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

use std::collections::HashMap;
use std::mem;

use log::*;

use crate::hg_utils::*;
use crate::hg_formats::*;
use crate::hg_math::*;


/* Various error messages */
const INVALID_K_BAR_MSG : &str = "Average degree must be greater than 0 and less than n-1.";
const INVALID_TEMPERATURE : &str = "Temperature must be positive (t >= 0).";
const INVALID_GAMMA : &str = "Gamma must be greater or equal 2 (gamma >= 2).";
const INVALID_GAMMA_WITH_ZETA : &str = "Zeta or eta make sense only at finite values of gamma.";


/* ================= graph construction utilities ================= */


/* RPrecomputedsinhcosh is a structure that contains the precomputed values
 * of sinh(zeta * r) and cosh(zeta * r) for each value of r. It is a map:
 *
 *         r => (sinh(zeta*r), cosh(zeta*r))
 *
 * this structure does not need to be exported outside this coding unit
 */

type RPrecomputedsinhcosh = HashMap<u64, (f64, f64)>;

// cast to usize to allow f64 as hash map key
fn f2u(f: f64) -> u64 {
  unsafe { mem::transmute(f) }
}

fn hg_assign_coordinates(
        nodes: &mut Vec<HgCoordinateType>,
        params: &HgParametersType,
        in_par: &HgAlgorithmParametersType,
        mut r_psc: Option<&mut RPrecomputedsinhcosh>,
        rnd_01: &mut FnMut() -> f64) {
  debug!("\tAssigning coordinates");

  match params.gtype {
    HgGraphType::HyperbolicRgg
    | HgGraphType::HyperbolicStandard
    | HgGraphType::SoftConfigurationModel => {
      for _ in 0..params.expected_n {
        let zeta = params.zeta_eta;
        let r = hg_quasi_uniform_radial_coordinate(in_par.radius, in_par.alpha, rnd_01);
        if let Some(ref mut r_psc) = r_psc {
          r_psc.insert(f2u(r), ((zeta * r).sinh(), (zeta * r).cosh()));
        }
        nodes.push(HgCoordinateType {
          r: r,
          theta: hg_uniform_angular_coordinate(rnd_01)
        });
      }
    },
    HgGraphType::AngularRgg
    | HgGraphType::SoftRgg
    | HgGraphType::ErdosRenyi => {
      for _ in 0..params.expected_n {
        nodes.push(HgCoordinateType {
          r: in_par.radius, // HG_INF_RADIUS,
          theta: hg_uniform_angular_coordinate(rnd_01)
        });
      }
    }
  }
}

/* ================= useful mathematical functions  ================= */


fn hg_get_rr_from_numerical_integration(
        params: &HgParametersType,
        p: &HgAlgorithmParametersType) -> f64 {
  return hg_get_rr(params, p);
}

fn hg_get_lambda_from_gauss_hypergeometric_function(
        params: &HgParametersType,
        p: &HgAlgorithmParametersType) -> f64 {
  return hg_get_lambda(params, p);
}


/* ================= single model graph generators  ================= */


fn hg_hyperbolic_distance_hyperbolic_rgg_standard_(
        zeta_eta: f64,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType,
        r_psc: Option<&RPrecomputedsinhcosh>) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }

  // if the nodes have the same angular coordinates
  // then we return the euclidean distance
  if node1.theta == node2.theta {
    return (node1.r - node2.r).abs();
  }

  // equation 13
  let zeta = zeta_eta;
  let delta_theta = HG_PI - (HG_PI - (node1.theta - node2.theta).abs()).abs();
  let (part1, part2) = if let Some(r_psc) = r_psc {
    let n1 = r_psc.get(&f2u(node1.r)).unwrap();
    let n2 = r_psc.get(&f2u(node2.r)).unwrap();

    (n1.1 * n2.1,
     n1.0 * n2.0 * delta_theta.cos())
  } else {
    ((zeta * node1.r).cosh() * (zeta * node2.r).cosh(),
     (zeta * node1.r).sinh() * (zeta * node2.r).sinh() * delta_theta.cos())
  };
  return (part1 - part2).acosh() / zeta;
}

pub fn hg_hyperbolic_distance_hyperbolic_rgg_standard(
        zeta_eta: f64,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  hg_hyperbolic_distance_hyperbolic_rgg_standard_(zeta_eta, node1, node2, None)
}

fn hg_connection_probability_hyperbolic_rgg(
          params: &HgParametersType,
          p: &HgAlgorithmParametersType,
          node1: &HgCoordinateType,
          node2: &HgCoordinateType,
          r_psc: Option<&RPrecomputedsinhcosh>) -> f64 {
  // equation 32: Heaviside function
  if hg_hyperbolic_distance_hyperbolic_rgg_standard_(params.zeta_eta, node1, node2, r_psc) <= p.radius {
    return 1.0;
  } else {
    return 0.0;
  }
}

pub fn hg_hyperbolic_rgg(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        exp_gamma: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(INVALID_K_BAR_MSG);
  }

  if exp_gamma < 2.0 {
    return Err(INVALID_GAMMA);
  }

  if exp_gamma >= HG_INF_GAMMA {
    return Err(INVALID_GAMMA_WITH_ZETA);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Hyperbolic Random Geometric Graph");

  let params = HgParametersType::new(n, k_bar,
    exp_gamma, 0.0 /* t = 0 */,
    zeta, HgGraphType::HyperbolicRgg);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = HgAlgorithmParametersType::new();
  p.alpha = 0.5 * zeta * (exp_gamma - 1.0);
  p.radius = hg_get_rr_from_numerical_integration(&params, &p);

  let mut r_psc = RPrecomputedsinhcosh::new();
  hg_assign_coordinates(&mut nodes, &params, &p, Some(&mut r_psc), rnd_01);

  debug!("\tInternal parameters:");
  debug!("\t\tAlpha: {}", p.alpha);
  debug!("\t\tRadius: {}", p.radius);

  debug!("\tCreating links");

  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_hyperbolic_rgg(&params, &p, &nodes[id], &nodes[other_id], Some(&r_psc)) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_connection_probability_hyperbolic_standard(
        params: &HgParametersType,
        p: &HgAlgorithmParametersType,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType,
        r_psc: Option<&RPrecomputedsinhcosh>) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }

  // equation 12: Fermi-Dirac function
  let zeta = params.zeta_eta;
  let t =  params.temperature;
  let x = hg_hyperbolic_distance_hyperbolic_rgg_standard_(params.zeta_eta, node1, node2, r_psc);
  let exponent = (1.0 / t) * (zeta / 2.0)  * (x - p.radius);

  return 1.0 / (exponent.exp() + 1.0);
}

pub fn hg_hyperbolic_standard(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        exp_gamma: f64,
        temperature: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(INVALID_K_BAR_MSG);
  }

  if exp_gamma < 2.0 {
    return Err(INVALID_GAMMA);
  }

  if exp_gamma >= HG_INF_GAMMA {
    return Err(INVALID_GAMMA_WITH_ZETA);
  }

  if temperature < 0.0 {
    return Err(INVALID_TEMPERATURE);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Hyperbolic Standard Graph\n");

  let params = HgParametersType::new(n, k_bar, exp_gamma, temperature,
    zeta, HgGraphType::HyperbolicStandard);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = HgAlgorithmParametersType::new();

  // alpha calculation. different for cold and hot regimes
  p.alpha = if temperature <= 1.0 {
    0.5 * zeta * (exp_gamma - 1.0)
  } else {
    0.5 * (zeta / temperature) * (exp_gamma - 1.0)
  };

  p.radius = hg_get_rr_from_numerical_integration(&params, &p);

  debug!("\tInternal parameters:");
  debug!("\t\tAlpha: {}", p.alpha);
  debug!("\t\tRadius: {}", p.radius);

  let mut r_psc = RPrecomputedsinhcosh::new();
  hg_assign_coordinates(&mut nodes, &params, &p, Some(&mut r_psc), rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_hyperbolic_standard(&params, &p, &nodes[id], &nodes[other_id], Some(&r_psc)) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

pub fn hg_hyperbolic_distance_scm(
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // curvature is infinite, so 1/zeta goes to zero
  return node1.r + node2.r;
}

fn hg_connection_probability_scm(
        p: &HgAlgorithmParametersType,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  // equation (39)
  let x = hg_hyperbolic_distance_scm(node1, node2);
  let exponent = (p.eta / 2.0) * (x - p.radius);
  return 1.0 / (exponent.exp() + 1.0);
}

pub fn hg_soft_configuration_model(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        exp_gamma: f64,
        eta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(INVALID_K_BAR_MSG);
  }

  if exp_gamma < 2.0 {
    return Err(INVALID_GAMMA);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Soft Configuration Model Graph\n");

  let params = HgParametersType::new(n, k_bar, exp_gamma,
    HG_INF_TEMPERATURE /* t = inf */,
    eta, HgGraphType::SoftConfigurationModel);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  // zeta goes to infinity
  //graph.zeta_eta = numeric_limits<double>::max( );
  let mut p = HgAlgorithmParametersType::new();
  p.alpha = 0.5 * p.eta * (exp_gamma - 1.0);
  p.eta = params.zeta_eta;
  p.radius = hg_get_rr_from_numerical_integration(&params, &p);

  debug!("\t\talpha: {}", p.alpha);
  debug!("\t\teta: {}", p.eta);
  debug!("\t\tradius: {}", p.radius);

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_scm(&p, &nodes[id], &nodes[other_id]) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

pub fn hg_hyperbolic_distance_angular_soft_rgg(
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // delta theta
  return HG_PI - (HG_PI - (node1.theta - node2.theta).abs()).abs();
}

fn hg_connection_probability_angular_rgg(
        params: &HgParametersType,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  // equation 55: Heaviside function
  if hg_hyperbolic_distance_angular_soft_rgg(node1, node2)
     <= (HG_PI * params.expected_degree / (params.expected_n as f64)) {
    return 1.0;
  }
  return 0.0;
}

pub fn hg_angular_rgg(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(INVALID_K_BAR_MSG);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Angular Random Geometric Graph\n");

  let params = HgParametersType::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    0.0 /* t = 0 */,
    zeta, HgGraphType::AngularRgg);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = HgAlgorithmParametersType::new();
  p.radius = HG_INF_RADIUS;

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_angular_rgg(&params, &nodes[id], &nodes[other_id]) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_connection_probability_soft_rgg(
        params: &HgParametersType,
        p: &HgAlgorithmParametersType,
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  let x = hg_hyperbolic_distance_angular_soft_rgg(node1, node2);
  // equation 46
  let beta = 1.0 / params.temperature;
  return 1.0 / (1.0 + p.c * (x / HG_PI).powf(beta));
}

pub fn hg_soft_rgg(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        temperature: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(INVALID_K_BAR_MSG);
  }

  if temperature < 0.0 {
    return Err(INVALID_TEMPERATURE);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Soft Random Geometric Graph\n");

  let params = HgParametersType::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    temperature /* t = 0 */,
    zeta, HgGraphType::SoftRgg);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = HgAlgorithmParametersType::new();
  p.radius = HG_INF_RADIUS;
  p.c = hg_get_lambda_from_gauss_hypergeometric_function(&params, &p);

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_soft_rgg(&params, &p, &nodes[id], &nodes[other_id]) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

pub fn hg_hyperbolic_distance_er(
        node1: &HgCoordinateType,
        node2: &HgCoordinateType) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // there is no "real distance", indeed!
  return 1.0;
}

fn hg_connection_probability_er(
        params: &HgParametersType,
        _node1: &HgCoordinateType,
        _node2: &HgCoordinateType) -> f64 {
  // connection probability is given
  // by equation 61
  return 1.0 / (1.0 + (params.expected_n as f64) / (params.expected_degree as f64));
}

pub fn hg_erdos_renyi(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n as f64 - 1.0)) {
    return Err(INVALID_K_BAR_MSG);
  }

  let mut nodes = Vec::<HgCoordinateType>::with_capacity(n);
  let mut links = Vec::<HgConnectionType>::with_capacity(n);

  debug!("-> Erdos-Renyi Graph\n");

  let params = HgParametersType::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    HG_INF_TEMPERATURE /* t = inf */,
    zeta, HgGraphType::ErdosRenyi);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = HgAlgorithmParametersType::new();
  p.radius = HG_INF_RADIUS;

  debug!("\t\tradius: {} (INF)", HG_INF_RADIUS);
  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_er(&params, &nodes[id], &nodes[other_id]) {
        links.push(HgConnectionType {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}
