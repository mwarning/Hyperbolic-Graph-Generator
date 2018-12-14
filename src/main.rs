
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
  let seed = 1;           // seed for the random number generator

  fn dstr(b: bool) -> &'static str {
    if b { "(default)" } else { "" }
  }

  println!("Parameters:");
  println!(" Number of nodes [n]: {} {}", n, dstr(n == 1000));
  println!(" Expected average degree [k]: {} {}", k_bar, dstr(k_bar == 10.0));

  print!(" Expected power-law exponent [g]: ");
  if exp_gamma >= HG_INF_GAMMA {
    println!("INF {}", dstr(exp_gamma == 2.0));
  } else {
    println!("{} {}", exp_gamma, dstr(exp_gamma == 2.0));
  }

  if (exp_gamma < HG_INF_GAMMA) && (t >= HG_INF_TEMPERATURE) {
    println!(" Ratio zeta/T [eta]: {} {}", zeta_eta, dstr(zeta_eta == 1.0));
  } else {
    println!(" Square root of curvature [z]: {} {}", zeta_eta, dstr(zeta_eta == 1.0));
  }

  print!(" Temperature [t]: ");
  if t >= HG_INF_TEMPERATURE {
    println!("INF {}", dstr(t == 0.0));
  } else {
    println!("{} {}", t, dstr(t == 0.0));
  }

  println!("Seed [s]: {} {}", seed, dstr(seed == 1));

  match hg_graph_generator(n, k_bar, exp_gamma, t, zeta_eta, seed) {
    Ok((nodes, links)) => {
      println!("nodes:");
      for node in &nodes {
        println!(" ({:.2}, {:.2})", node.r, node.theta);
      }

      println!("links:");
      for link in &links {
        println!(" {} <-> {}", link.id, link.other_id);
      }

      println!("nodes: {}, links: {}", nodes.len(), links.len());
    },
    Err(msg) => {
      error!("Error: {}", msg);
    }
  }
}
