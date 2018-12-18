
#[macro_use]
mod hg_debug;
mod hg_formats;
mod hg_math;
mod hg_utils;
mod hg_gen_algorithms;
mod hg_graphs_lib;

use crate::hg_formats::*;
use crate::hg_graphs_lib::*;

use std::str::FromStr;
use std::env;


const NUM_DEFAULT : usize = 100;
const EXP_DEGREE_DEFAULT : f64 = 10.0;
const EXP_GAMMA_DEFAULT : f64 = 2.0;
const TEMP_DEFAULT : f64 = 0.0;
const ZETA_ETA_DEFAULT : f64 = 1.0;
const SEED_DEFAULT : u32 = 1;
const OUTPUT_FORMAT_DEFAULT : &'static str = "tsv";


fn usage() {
  println!(
    concat!("Generate different types of undirected hyperbolic graphs.\n",
    "\n",
    "PARAMETERS:\n",
    "\t-n\tNumber of nodes to generate (default: {})\n",
    "\t-k\tExpected average degree (default: {})\n",
    "\t-g\tExpected gamma or gamma out (default: {})\n",
    "\t-t\tTemperature (default: {})\n",
    "\t-z\tParameter associated with curvature (default: {})\n",
    "\t-s\tRandom generator seed (default: {})\n",
    "\t-f\tOutput format tsv/json (default: {})\n",
    "\t-h\tPrint this help\n",
    "\n",
    "DESCRIPTION:\n",
    "\tAdapted description of the C++ version:\n",
    "\tThe program generates a graph that describes the geometric coordinates and links\n",
    "\tof a hyperbolic graph compatible with the parameters provided by the user.\n",
    "\tThe program generates random hyperbolic graphs according to the models in:\n",
    "\t  * http://dx.doi.org/10.1103/PhysRevE.82.036106\n",
    "\tA description of how the hyperbolic graph generator works can be found at:\n",
    "\t  * http://arxiv.org/abs/1503.05180\n",
    "\n",
    "OUTPUT:\n",
    "\tThe program prints TSV or JSON output to the standard output.\n",
    "\tTSV output is structured the folowing way:\n",
    "\n",
    "\t# The first line describes the main graph parameters:\n",
    "\tN <num nodes> T <temperature> G <gamma> K <avg.degree> Z <zeta> S <seed> I <initial_node_id>\n",
    "\n",
    "\t# A line for each node and its polar coordinates:\n",
    "\t<node_id>\t<radial coordinate>\t<angular coordinate>\n",
    "\n",
    "\t# A line for each link:\n",
    "\t<node_id>\t<node id>\n"
    ),
      NUM_DEFAULT, EXP_DEGREE_DEFAULT, EXP_GAMMA_DEFAULT,
      TEMP_DEFAULT, ZETA_ETA_DEFAULT, SEED_DEFAULT, OUTPUT_FORMAT_DEFAULT
    );
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
  let mut num = NUM_DEFAULT;               // number of nodes in the graph
  let mut exp_degree = EXP_DEGREE_DEFAULT; // expected average degree
  let mut exp_gamma = EXP_GAMMA_DEFAULT;   // expected gamma or gamma out
  let mut temp = TEMP_DEFAULT;             // temperature
  let mut zeta_eta = ZETA_ETA_DEFAULT;     // parameter associated with curvature
  let mut seed = SEED_DEFAULT;             // seed for the random number generator
  let mut output_format = OUTPUT_FORMAT_DEFAULT.to_string();

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
      }
    }
    i += 1;
  }

  // get default string
  fn dstr(b: bool) -> &'static str {
    if b { "\t(default)" } else { "" }
  }

  // get comma
  fn dcom(b: bool) -> &'static str {
    if b { "," } else { "" }
  }

  eprintln!("Parameters:");
  eprintln!("\tNumber of nodes [n]:\t\t\t{}{}", num, dstr(num == 1000));
  eprintln!("\tExpected average degree [k]:\t\t{}{}", exp_degree, dstr(exp_degree == 10.0));

  eprint!("\tExpected power-law exponent [g]:\t");
  if exp_gamma >= HG_INF_GAMMA {
    eprintln!("INF {}", dstr(exp_gamma == 2.0));
  } else {
    eprintln!("{} {}", exp_gamma, dstr(exp_gamma == 2.0));
  }

  if (exp_gamma < HG_INF_GAMMA) && (temp >= HG_INF_TEMPERATURE) {
    eprintln!("\tRatio zeta/T [eta]:\t\t{}{}", zeta_eta, dstr(zeta_eta == 1.0));
  } else {
    eprintln!("\tSquare root of curvature [z]:\t\t{}{}", zeta_eta, dstr(zeta_eta == 1.0));
  }

  eprint!("\tTemperature [t]:\t\t\t");
  if temp >= HG_INF_TEMPERATURE {
    eprintln!("INF{}", dstr(temp == 0.0));
  } else {
    eprintln!("{}{}", temp, dstr(temp == 0.0));
  }

  eprintln!("\tSeed [s]:\t\t\t\t{}{}\n", seed, dstr(seed == 1));

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
          println!("  {{\"x\": {:.10}, \"y\": {:.10}}}{}", node.r, node.theta, dcom((i + 1) != nlen));
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
        let starting_id = 1;
        print!("N\t{}\tT\t{:.10}\tG\t{:.10}\tK\t{:.10}\t", num, temp, exp_gamma, exp_degree);
        if (temp >= HG_INF_TEMPERATURE) && (exp_gamma < HG_INF_GAMMA) {
          print!("eta\t{:.10}\t", zeta_eta);
        } else {
          print!("Z\t{:.10}\t", zeta_eta);
        }
        println!("S\t{}\tI {}", seed, starting_id);

        for (i, node) in nodes.iter().enumerate() {
          println!("{}\t{:.10}\t{:.10}", i + starting_id, node.r, node.theta);
        }

        for link in &links {
          println!("{}\t{}", link.id + starting_id, link.other_id + starting_id);
        }
      }
      eprintln!("Generated:\t{} links", links.len());
    },
    Err(msg) => {
      eprintln!("Error: {}", msg);
    }
  }
}
