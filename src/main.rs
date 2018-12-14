
mod hg_formats;
mod hg_math;
mod hg_utils;
mod hg_gen_algorithms;
mod hg_graphs_lib;

use crate::hg_formats::*;
use crate::hg_graphs_lib::*;

use log::*;


fn main() {
  /* HG Graph Generator: default parameters */
  let n = 1000;           // number of nodes in the graph
  let k_bar = 10.0f64;    // expected average degree
  let exp_gamma = 2.0f64; // expected gamma or gamma out
  let t = 0.0f64;         // temperature
  let zeta_eta = 1.0f64;  // parameter associated with curvature

  let default = "(default)";
  let no_default = "";

  let n_default_str = if n == 1000 { default } else { no_default };
  let k_bar_default_str = if k_bar == 10.0 { default } else { no_default };
  let exp_gamma_default_str = if exp_gamma == 2.0 { default } else { no_default };
  let zeta_eta_default_str = if zeta_eta == 1.0 { default } else { no_default };
  let t_default_str = if t == 0.0 { default } else { no_default };

  println!("Parameters:");
  println!(" Number of nodes [n]: {} {}", n, n_default_str);
  println!(" Expected average degree [k]: {} {}", k_bar, k_bar_default_str);

  print!(" Expected power-law exponent [g]: ");
  if exp_gamma >= HG_INF_GAMMA {
    print!("INF {}", exp_gamma_default_str);
  } else {
    print!("{} {}", exp_gamma, exp_gamma_default_str);
  }

  if (exp_gamma < HG_INF_GAMMA) && (t >= HG_INF_TEMPERATURE) {
    println!(" Ratio zeta/T [eta]: {} {}", zeta_eta, zeta_eta_default_str);
  } else {
    println!(" Square root of curvature [z]: {} {}", zeta_eta, zeta_eta_default_str);
  }

  print!(" Temperature [t]: ");
  if (t >= HG_INF_TEMPERATURE) {
    println!("INF {}", t_default_str);
  } else {
    println!("{} {}", t, t_default_str);
  }

  match hg_graph_generator(n, k_bar, exp_gamma, t, zeta_eta) {
    Ok((nodes, links)) => {
      // do stuff
    },
    Err(msg) => {
      error!("Error: {}", msg);
    }
  }
}
