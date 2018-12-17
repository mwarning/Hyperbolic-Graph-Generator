
mod hg_formats;
mod hg_math;
mod hg_utils;
mod hg_gen_algorithms;
mod hg_graphs_lib;

use crate::hg_formats::*;
use crate::hg_graphs_lib::*;

use std::str::FromStr;
use std::env;
use log::*;


fn usage() {
  println!(
    concat!("Generate different types of undirected hyperbolic graphs.\n",
    "\n",
    " -n\tNumber of nodes to generate\n",
    " -k\tExpected average degree\n",
    " -g\tExpected gamma or gamma out\n",
    " -t\tTemperature\n",
    " -z\tParameter associated with curvature\n",
    " -s\tRandom generator seed\n",
    " -f\tOutput format (tsv or json)\n",
    " -h\tPrint this help\n",
    "\n",
    " Adapted description of the C++ version:\n",
    "  The program generates a graph that describes the geometric coordinates and links\n",
    "  of a hyperbolic graph compatible with the parameters provided by the user.\n",
    "  The program generates random hyperbolic graphs according to the models in:\n",
    "    * http://dx.doi.org/10.1103/PhysRevE.82.036106\n",
    "  A description of how the hyperbolic graph generator works can be found at:\n",
    "    * http://arxiv.org/abs/1503.05180\n"
    ));
}

fn get_arg<T : FromStr>(i: usize, args: &Vec<String>) -> T {
  let name = &args[i];

  if (i+1) >= args.len() {
    eprintln!("Missing argument for {}", name);
    std::process::exit(1);
  }

  let value = &args[i+1];
  if let Ok(val) = value.parse::<T>() {
    return val;
  } else {
    eprintln!("Invalid value for {}: {}", name, value);
    std::process::exit(1);
  }
}

fn main() {
  /* HG Graph Generator: default parameters */
  let mut num = 100;            // number of nodes in the graph
  let mut exp_degree = 10.0f64; // expected average degree
  let mut exp_gamma = 2.0f64;   // expected gamma or gamma out
  let mut temp = 0.0f64;        // temperature
  let mut zeta_eta = 1.0f64;    // parameter associated with curvature
  let mut seed = 1;             // seed for the random number generator
  let mut output_format = "tsv".to_string();

  let args = env::args().collect::<Vec<String>>();
  let mut i = 1;
  while i < args.len() {
    match args[i].as_ref() {
      "-n" => {
        num = get_arg(i, &args);
        i += 1;
      },
      "-k" => {
        exp_degree = get_arg(i, &args);
        i += 1;
      },
      "-g" => {
        exp_gamma = get_arg(i, &args);
        i += 1;
      },
      "-t" => {
        temp = get_arg(i, &args);
        i += 1;
      },
      "-z" => {
        zeta_eta = get_arg(i, &args);
        i += 1;
        //zeta_eta_provided = true; //needed?
      },
      "-s" => {
        seed = get_arg(i, &args);
        i += 1;
      },
      "-f" => {
        output_format = get_arg(i, &args);
        i += 1;
      },
      "-h" => {
        usage();
        return;
      },
      _ => {
        eprintln!("Unknown argument: {}", args[i]);
        std::process::exit(1);
        //Err((), )
      }
    }
    i += 1;
  }

  // get default string
  fn dstr(b: bool) -> &'static str {
    if b { "(default)" } else { "" }
  }

  // get comma
  fn dcom(b: bool) -> &'static str {
    if b { "," } else { "" }
  }

  eprintln!("Parameters:");
  eprintln!(" Number of nodes [n]: {} {}", num, dstr(num == 1000));
  eprintln!(" Expected average degree [k]: {} {}", exp_degree, dstr(exp_degree == 10.0));

  eprint!(" Expected power-law exponent [g]: ");
  if exp_gamma >= HG_INF_GAMMA {
    eprintln!("INF {}", dstr(exp_gamma == 2.0));
  } else {
    eprintln!("{} {}", exp_gamma, dstr(exp_gamma == 2.0));
  }

  if (exp_gamma < HG_INF_GAMMA) && (temp >= HG_INF_TEMPERATURE) {
    eprintln!(" Ratio zeta/T [eta]: {} {}", zeta_eta, dstr(zeta_eta == 1.0));
  } else {
    eprintln!(" Square root of curvature [z]: {} {}", zeta_eta, dstr(zeta_eta == 1.0));
  }

  eprint!(" Temperature [t]: ");
  if temp >= HG_INF_TEMPERATURE {
    eprintln!("INF {}", dstr(temp == 0.0));
  } else {
    eprintln!("{} {}", temp, dstr(temp == 0.0));
  }

  eprintln!(" Seed [s]: {} {}", seed, dstr(seed == 1));

  match hg_graph_generator(num, exp_degree, exp_gamma, temp, zeta_eta, seed) {
    Ok((nodes, links)) => {
      if output_format == "json" {
        // JSON output
        println!("{{");

        println!(" \"parameters\": {{");
        println!("  \"exp_degree\": {},", exp_degree);
        println!("  \"exp_gamma\": {},", exp_gamma);
        println!("  \"temperature\": {},", temp);
        println!("  \"zeta_eta\": {},", zeta_eta);
        println!("  \"seed\": {}", seed);
        println!(" }},");

        let nlen = nodes.len();
        println!(" \"nodes\": [");
        for (i, node) in nodes.iter().enumerate() {
          println!("  {{\"x\": {}, \"y\": {}}}{}", node.r, node.theta, dcom((i + 1) != nlen));
        }
        println!(" ],");

        let llen = links.len();
        println!(" \"links\": [");
        for (i, link) in links.iter().enumerate() {
          println!("  {{\"source\": {}, \"target\": {}}}{}", link.id, link.other_id, dcom((i + 1) != llen));
        }
        println!(" ]");

        println!("}}");
      } else {
        // TSV output
        print!("#N\t{}\tT\t{}\tG\t{}\tK\t{}\t", num, temp, exp_gamma, exp_degree);
        if (temp >= HG_INF_TEMPERATURE) && (exp_gamma < HG_INF_GAMMA) {
          print!("eta\t{}\t", zeta_eta);
        } else {
          print!("Z\t{}\t", zeta_eta);
        }
        println!("S\t{}", seed);

        for node in &nodes {
          println!("{}\t{}", node.r, node.theta);
        }

        for link in &links {
          println!("{}\t{}", link.id, link.other_id);
        }
      }
    },
    Err(msg) => {
      error!("Error: {}", msg);
    }
  }
}
