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
const invalid_k_bar_msg : &str = "Average degree must be greater than 0 and less than n-1.";
const invalid_temperature : &str = "Temperature must be positive (t >= 0).";
const invalid_gamma : &str = "Gamma must be greater or equal 2 (gamma >= 2).";
const invalid_gamma_with_zeta : &str = "Zeta or eta make sense only at finite values of gamma.";


/* ================= graph construction utilities ================= */


/* r_precomputedsinhcosh is a structure that contains the precomputed values
 * of sinh(zeta * r) and cosh(zeta * r) for each value of r. It is a map:
 *
 *         r => (sinh(zeta*r), cosh(zeta*r))
 *
 * this structure does not need to be exported outside this coding unit
 */

type r_precomputedsinhcosh = HashMap<u64, (f64, f64)>;

// allow f64 as hash map key
fn f2u(f: f64) -> u64 {
  unsafe { mem::transmute(f) }
}

fn hg_assign_coordinates(
        nodes: &mut Vec<hg_coordinate_t>,
        params: &hg_parameters_t,
        in_par: &hg_algorithm_parameters_t,
        mut r_psc: Option<&mut r_precomputedsinhcosh>,
        rnd_01: &mut FnMut() -> f64) {
  debug!("\tAssigning coordinates");

  match params.gtype {
    hg_graph_type::HYPERBOLIC_RGG
    | hg_graph_type::HYPERBOLIC_STANDARD
    | hg_graph_type::SOFT_CONFIGURATION_MODEL => {
      for _ in 0..params.expected_n {
        let zeta = params.zeta_eta;
        let r = hg_quasi_uniform_radial_coordinate(in_par.radius, in_par.alpha, rnd_01);
        if let Some(ref mut r_psc) = r_psc {
          r_psc.insert(f2u(r), ((zeta * r).sinh(), (zeta * r).cosh()));
        }
        nodes.push(hg_coordinate_t {
          r: r,
          theta: hg_uniform_angular_coordinate(rnd_01)
        });
      }
    },
    hg_graph_type::ANGULAR_RGG
    | hg_graph_type::SOFT_RGG
    | hg_graph_type::ERDOS_RENYI => {
      for _ in 0..params.expected_n {
        nodes.push(hg_coordinate_t {
          r: in_par.radius, // HG_INF_RADIUS,
          theta: hg_uniform_angular_coordinate(rnd_01)
        });
      }
    }
  }
}

/*
fn hg_init_graph(n: u32,
                 k_bar: f64,
                 exp_gamma: f64,
                 t: f64,
                 zeta_eta: f64,
                 seed: i32,
                 gt: hg_graph_type) -> hg_graph_t {
  debug!("\tGraph initialization");
  /* initialize the graph structure with the 
   * parameters provided in input by the user */

  // Init random generator
  HG_Random::init(seed); // TODO: move out and create hg_graph_t::new()

  hg_graph_t {
    gtype: gt,
    expected_n: n,
    temperature: t,
    expected_gamma: exp_gamma,
    expected_degree: k_bar,
    zeta_eta: zeta_eta,
    seed: seed,
    starting_id: 1,
    data: vec![],
  }
}*/


/* ================= useful mathematical functions  ================= */


fn hg_get_R_from_numerical_integration(
        params: &hg_parameters_t,
        p: &hg_algorithm_parameters_t) -> f64 {
  return hg_get_R(params, p);
}

fn hg_get_lambda_from_Gauss_hypergeometric_function(
        params: &hg_parameters_t,
        p: &hg_algorithm_parameters_t) -> f64 {
  return hg_get_lambda(params, p);
}


/* ================= single model graph generators  ================= */


fn hg_hyperbolic_distance_hyperbolic_rgg_standard(
        zeta_eta: f64,
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t,
        r_psc: Option<&r_precomputedsinhcosh>) -> f64 {
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
    println!("hit cache!");
    let n1 = r_psc.get(&f2u(node1.r)).unwrap();
    let n2 = r_psc.get(&f2u(node2.r)).unwrap();

    (n1.1 * n2.1,
     n1.0 * n2.0 * delta_theta.cos())

    //(r_psc[node1.r].1 * r_psc[node2.r].1,
    // r_psc[node1.r].0 * r_psc[node2.r].0 * delta_theta.cos())
  } else {    
    ((zeta * node1.r).cosh() * (zeta * node2.r).cosh(),
     (zeta * node1.r).sinh() * (zeta * node2.r).sinh() * delta_theta.cos())
  };
  return (part1 - part2).acosh() / zeta;
}

fn hg_hyperbolic_distance_hyperbolic_rgg_std(
        zeta_eta: f64,
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
  hg_hyperbolic_distance_hyperbolic_rgg_standard(zeta_eta, node1, node2, None)
}

fn hg_connection_probability_hyperbolic_rgg(
          params: &hg_parameters_t,
          p: &hg_algorithm_parameters_t,
          node1: &hg_coordinate_t,
          node2: &hg_coordinate_t,
          r_psc: Option<&r_precomputedsinhcosh>) -> f64 {
  // equation 32: Heaviside function
  if hg_hyperbolic_distance_hyperbolic_rgg_standard(params.zeta_eta, node1, node2, r_psc) <= p.radius {
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
    return Err(invalid_k_bar_msg);
  }

  if exp_gamma < 2.0 {
    return Err(invalid_gamma);
  }

  if exp_gamma >= HG_INF_GAMMA {
    return Err(invalid_gamma_with_zeta);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Hyperbolic Random Geometric Graph");

  let params = hg_parameters_t::new(n, k_bar,
    exp_gamma, 0.0 /* t = 0 */,
    zeta, hg_graph_type::HYPERBOLIC_RGG);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = hg_algorithm_parameters_t::new();
  p.alpha = 0.5 * zeta * (exp_gamma - 1.0);
  p.radius = hg_get_R_from_numerical_integration(&params, &p);

  let mut r_psc = r_precomputedsinhcosh::new(); 
  hg_assign_coordinates(&mut nodes, &params, &p, Some(&mut r_psc), rnd_01);

  debug!("\tInternal parameters:");
  debug!("\t\tAlpha: {}", p.alpha);
  debug!("\t\tRadius: {}", p.radius);

  debug!("\tCreating links");

  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_hyperbolic_rgg(&params, &p, &nodes[id], &nodes[other_id], Some(&r_psc)) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_connection_probability_hyperbolic_standard(
        params: &hg_parameters_t,
        p: &hg_algorithm_parameters_t,
        node1: &hg_coordinate_t,
        node2: &hg_coordinate_t,
        r_psc: Option<&r_precomputedsinhcosh>) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }

  // equation 12: Fermi-Dirac function
  let zeta = params.zeta_eta;
  let t =  params.temperature;
  let x = hg_hyperbolic_distance_hyperbolic_rgg_standard(params.zeta_eta, node1, node2, r_psc);
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
    return Err(invalid_k_bar_msg);
  }

  if exp_gamma < 2.0 {
    return Err(invalid_gamma);
  }

  if exp_gamma >= HG_INF_GAMMA {
    return Err(invalid_gamma_with_zeta);
  }

  if temperature < 0.0 {
    return Err(invalid_temperature);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Hyperbolic Standard Graph\n");

  let params = hg_parameters_t::new(n, k_bar, exp_gamma, temperature,
    zeta, hg_graph_type::HYPERBOLIC_STANDARD);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = hg_algorithm_parameters_t::new();

  // alpha calculation. different for cold and hot regimes
  p.alpha = if temperature <= 1.0 {
    0.5 * zeta * (exp_gamma - 1.0)
  } else {
    0.5 * (zeta / temperature) * (exp_gamma - 1.0)
  };

  p.radius = hg_get_R_from_numerical_integration(&params, &p);

  debug!("\tInternal parameters:");
  debug!("\t\tAlpha: {}", p.alpha);
  debug!("\t\tRadius: {}", p.radius);

  let mut r_psc = r_precomputedsinhcosh::new(); 
  hg_assign_coordinates(&mut nodes, &params, &p, Some(&mut r_psc), rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_hyperbolic_standard(&params, &p, &nodes[id], &nodes[other_id], Some(&r_psc)) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_hyperbolic_distance_scm(
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // curvature is infinite, so 1/zeta goes to zero  
  return node1.r + node2.r;
}

fn hg_connection_probability_scm(
        p: &hg_algorithm_parameters_t,
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
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
    return Err(invalid_k_bar_msg);
  }

  if exp_gamma < 2.0 {
    return Err(invalid_gamma);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Soft Configuration Model Graph\n");

  let params = hg_parameters_t::new(n, k_bar, exp_gamma,
    HG_INF_TEMPERATURE /* t = inf */,
    eta, hg_graph_type::SOFT_CONFIGURATION_MODEL);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  // zeta goes to infinity
  //graph.zeta_eta = numeric_limits<double>::max( );
  let mut p = hg_algorithm_parameters_t::new();
  p.alpha = 0.5 * p.eta * (exp_gamma - 1.0);
  p.eta = params.zeta_eta;
  p.radius = hg_get_R_from_numerical_integration(&params, &p);

  debug!("\t\talpha: {}", p.alpha);
  debug!("\t\teta: {}", p.eta);
  debug!("\t\tradius: {}", p.radius);

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_scm(&p, &nodes[id], &nodes[other_id]) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_hyperbolic_distance_angular_soft_rgg(
        node1: &hg_coordinate_t,
        node2: &hg_coordinate_t) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // delta theta
  return HG_PI - (HG_PI - (node1.theta - node2.theta).abs()).abs();
}

fn hg_connection_probability_angular_rgg(
        params: &hg_parameters_t,
        node1: &hg_coordinate_t,
        node2: &hg_coordinate_t) -> f64 {
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
    return Err(invalid_k_bar_msg);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Angular Random Geometric Graph\n");

  let params = hg_parameters_t::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    0.0 /* t = 0 */,
    zeta, hg_graph_type::ANGULAR_RGG);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = hg_algorithm_parameters_t::new();
  p.radius = HG_INF_RADIUS;

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_angular_rgg(&params, &nodes[id], &nodes[other_id]) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_connection_probability_soft_rgg(
        params: &hg_parameters_t,
        p: &hg_algorithm_parameters_t,
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
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
    return Err(invalid_k_bar_msg);
  }

  if temperature < 0.0 {
    return Err(invalid_temperature);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Soft Random Geometric Graph\n");

  let params = hg_parameters_t::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    temperature /* t = 0 */,
    zeta, hg_graph_type::SOFT_RGG);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = hg_algorithm_parameters_t::new();
  p.radius = HG_INF_RADIUS;
  p.c = hg_get_lambda_from_Gauss_hypergeometric_function(&params, &p);

  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {
      if rnd_01() < hg_connection_probability_soft_rgg(&params, &p, &nodes[id], &nodes[other_id]) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}

fn hg_hyperbolic_distance_er(
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
  // check if it is the same node
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }
  // there is no "real distance", indeed!
  return 1.0;
}

fn hg_connection_probability_er(
        params: &hg_parameters_t,
        node1: &hg_coordinate_t, 
        node2: &hg_coordinate_t) -> f64 {
  // connection probability is given
  // by equation 61
  return 1.0 / (1.0 + (params.expected_n as f64) / (params.expected_degree as f64)); 
}

pub fn hg_erdos_renyi(
        n: usize,
        rnd_01: &mut FnMut() -> f64,
        k_bar: f64,
        zeta: f64) -> Result<Graph, &'static str>{

  if (k_bar < 1.0) || (k_bar > (n - 1) as f64) {
    return Err(invalid_k_bar_msg);
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();

  debug!("-> Erdos-Renyi Graph\n");

  let params = hg_parameters_t::new(n, k_bar,
    HG_INF_GAMMA /* exp_gamma = inf */,
    HG_INF_TEMPERATURE /* t = inf */,
    zeta, hg_graph_type::ERDOS_RENYI);

  // computing internal parameters
  debug!("\tInternal parameters computation");
  let mut p = hg_algorithm_parameters_t::new();
  p.radius = HG_INF_RADIUS;

  debug!("\t\tradius: {} (INF)", HG_INF_RADIUS);
  hg_assign_coordinates(&mut nodes, &params, &p, None, rnd_01);

  debug!("\tCreating links");
  for id in 0..params.expected_n {
    for other_id in (id + 1)..params.expected_n {  
      if rnd_01() < hg_connection_probability_er(&params, &nodes[id], &nodes[other_id]) {
        links.push(hg_connection_t {id: id, other_id: other_id});
      }
    }
  }

  Ok((nodes, links))
}


/* ================= hyperbolic distance function  ================= */


fn hg_hyperbolic_distance(
        params: &hg_parameters_t,
        node1: &hg_coordinate_t,
        node2: &hg_coordinate_t) -> f64 {
  if (node1.r == node2.r) && (node1.theta == node2.theta) {
    return 0.0;
  }

  match params.gtype {
    hg_graph_type::HYPERBOLIC_RGG | hg_graph_type::HYPERBOLIC_STANDARD => 
      hg_hyperbolic_distance_hyperbolic_rgg_std(params.zeta_eta, node1, node2),
    hg_graph_type::SOFT_CONFIGURATION_MODEL =>
      hg_hyperbolic_distance_scm(node1, node2),
    hg_graph_type::ANGULAR_RGG | hg_graph_type::SOFT_RGG =>
      hg_hyperbolic_distance_angular_soft_rgg(node1, node2),
    hg_graph_type::ERDOS_RENYI =>
      hg_hyperbolic_distance_er(node1, node2)
  }
}
