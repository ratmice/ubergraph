## Ubergraph
A small prototype based on

- [Ubergraphs](https://arxiv.org/pdf/1704.05547.pdf) a definition of recursive hypergraph structure


Adjacency representation is currently done via a vector of sets.
No algorithms have been implemented yet.

- [X] Adding edges and nodes
- [X] Conversion of Levi or Incidence graph into petgraph.
- [ ] Everything else

## Notes, and initial impressions of the implementation

There are lots of rough edges (ownership, Node and edge Index types),
in particular due to rusts lack of tail recursion unrolling some of these instances of recursion
to avoid stack frames on recursive calls seems like it is going to be painful.
