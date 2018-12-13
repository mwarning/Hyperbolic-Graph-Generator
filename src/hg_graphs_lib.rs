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

/*
hg_graph_t * hg_read_graph(const string filename) {
  // graph_t pointer
  hg_graph_t *g = NULL;
  // graph generation parameters
  string par; // name of parameter
  // node id
  int i;
  int node1,node2;
  // coordinates
  double radial;
  double angular;
  // opening file
  ifstream file;
  int expected_n;

  file.open(filename.c_str(), ios::in);
  if(!file.fail() and file.is_open()) {
    // reading first line
    file >> par >> expected_n;
    g = new hg_graph_t(expected_n);
    g.expected_n = expected_n;
    file >> par >> g.temperature;
    file >> par >> g.expected_gamma;
    file >> par >> g.expected_degree;
    file >> par >> g.zeta_eta;
    file >> par >> g.seed;
    file >> par >> g.starting_id;
    // reading coordinates
    for(i = 0; i < expected_n; i++) {
      file >> node1 >> radial >> angular;
      (*g)[i].r = radial;
      (*g)[i].theta = angular;
    }
    // reading links
    unsigned int starting_id = g.starting_id;
    while(file >> node1 >> node2) {
      add_edge(node1-starting_id, node2-starting_id, *g);
    }
  }
  else {
    // warning
    hg_log_err("File %s cannot be opened", filename.c_str());
    return NULL;
  }
  file.close();
  // infer graph type from parameters 
  g.type = hg_infer_hg_type(g);
  return g; 
}

void hg_print_graph(const hg_graph_t *g, const string filename) {

  if(g == NULL) {
    hg_log_err("Warning: empty data structure, no file written");
    return;
  }
  ofstream file;

  file.open(filename.c_str(), ios::out);
  if(!file.fail() and file.is_open()) {
    file << std::setprecision(10) << std::fixed;
    // hg_graph parameters
    file << "N" << "\t" << g.expected_n << "\t";
    file << "T" << "\t" << g.temperature << "\t";
    file << "G" << "\t" << g.expected_gamma << "\t";
    file << "K" << "\t" << g.expected_degree << "\t";
    if(g.temperature >= HG_INF_TEMPERATURE &&
       g.expected_gamma < HG_INF_GAMMA) {
      file << "eta" << "\t" << g.zeta_eta << "\t";
    }
    else {
      file << "Z" << "\t" << g.zeta_eta << "\t";
    }
    file << "S" << "\t" << g.seed << "\t";
    file << "I" << "\t" << g.starting_id << endl;
    // hg_graph vertex coordinates
    unsigned int starting_id = g.starting_id;
    hg_graph_t::vertex_iterator vertexIt, vertexEnd;
    boost::tie(vertexIt, vertexEnd) = vertices(*g);
    for (; vertexIt != vertexEnd; ++vertexIt) { 
      file << *vertexIt + starting_id << "\t";
      file << (*g)[*vertexIt].r << "\t";
      file << (*g)[*vertexIt].theta << endl;
    }
    // hg_graph edgelist
    hg_graph_t::edge_iterator edgeIt, edgeEnd;
    boost::tie(edgeIt, edgeEnd) = edges(*g);
    for (; edgeIt != edgeEnd; ++edgeIt) { 
      file << source(*edgeIt, *g) + starting_id << "\t";
      file << target(*edgeIt, *g) + starting_id << endl;
    }
  }
  file.close();
  return;
}

void hg_init_random_generator(const unsigned int seed) {
  HG_Random::init(seed);
}


double hg_rand_01_wrapper() {
  // return =  ((double) rand() / (RAND_MAX));
  return  HG_Random::get_random_01_value();
}
*/

pub fn hg_graph_generator(
        n: usize,
        k_bar: f64, 
				exp_gamma: f64,
        t : f64, 
				zeta: f64) -> Graph {
  let mut rnd_01 = || { 0.0f64 };//and set seed
  let gt = hg_infer_hg_type(exp_gamma, t);

  /*
   //Limit cases that we do not take into account
  if(n < 3){
    hg_enduser_warning("Number of nodes must be n>=3. \n\t  Quitting.");
    return 1;
  }
  if (k_bar < 1) || (k_bar > n-1){
    hg_enduser_warning("Avg degree must be greater than 0 and less than n-1. \n\t  Quitting.");   
    return 1;
  }
  if t < 0 {
    hg_enduser_warning("Temperature must be positive (t >= 0). \n\t  Quitting.");
    return 1;
  }  
  if exp_gamma < 2.0 {
    hg_enduser_warning("Gamma must be greater or equal 2 (Gamma >= 2). \n\t  Quitting.");
    return 1;
  }  

  // Warnings 
  //if(zeta_provided and (t>=HG_INF_TEMPERATURE || exp_gamma>=HG_INF_GAMMA)) {
  if (zeta_eta_provided && exp_gamma>=HG_INF_GAMMA) {
    hg_enduser_warning("zeta or eta make sense only at finite values of gamma. \n\t  The provided value of zeta (or eta) will be ignored.");
    zeta_eta = 1;
  }*/

  if let Ok(graph) = match gt {
    HYPERBOLIC_RGG =>
      hg_hyperbolic_rgg(n, &mut rnd_01, k_bar, exp_gamma, zeta),
    HYPERBOLIC_STANDARD => 
      hg_hyperbolic_standard(n, &mut rnd_01, k_bar, exp_gamma, t, zeta),
    SOFT_CONFIGURATION_MODEL =>
      hg_soft_configuration_model(n, &mut rnd_01, k_bar, exp_gamma, zeta),
    ANGULAR_RGG =>
      hg_angular_rgg(n, &mut rnd_01, k_bar, zeta),
    SOFT_RGG =>
      hg_soft_rgg(n, &mut rnd_01, k_bar, t, zeta),
    ERDOS_RENYI =>
      hg_erdos_renyi(n, &mut rnd_01, k_bar, zeta),
  } {
    
  }

  let mut nodes = Vec::<hg_coordinate_t>::new();
  let mut links = Vec::<hg_connection_t>::new();
  (nodes, links)
}
