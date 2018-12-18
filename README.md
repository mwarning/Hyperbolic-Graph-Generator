# Hyperbolic-Graph-Generator


This is a Rust port of the [Hyperbolic-Graph-Generator 1.0.3](https://github.com/named-data/Hyperbolic-Graph-Generator) excluding most of the tools.

The program generates a graph that describes the geometric coordinates and the links of a hyperbolic graph compatible with the parameters provided by the user. The program generates random hyperbolic graphs according to the models in: http://dx.doi.org/10.1103/PhysRevE.82.036106

A description of how the hyperbolic graph generator works can be found at: http://arxiv.org/abs/1503.05180

TODO:
- crate polish
- better API
- c bindings

## Differences to C++ Version

- twice as fast
- the graph_properties and greedy_routing tool was not ported
- slightly different API
  - C++ uses boost::adjacency_list for the graph representation. In Rust we use two Vectors of node and link objects.
- optional json output and the original tsv output
