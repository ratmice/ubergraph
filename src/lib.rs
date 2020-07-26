use either::Either;
use petgraph::graph::Graph;
use petgraph::graph::NodeIndex;
use petgraph::Directed;
use petgraph::Direction;

use std::iter::Iterator;
use std::marker::PhantomData;

///! A Simple stupid attempt at building ubergraphs based on vectors.
///! While reusing a bunch of types from petgraph.
///!
///! Probably not ideal, as it means you cannot remove nodes edges, or edge members without
///! affecting the stability of indices.
///!
///! Just supports the basics of adding vertices and edges, and computing the incidence graph.

#[derive(Copy, Clone)]
pub enum EdgeMember<VIx, EIx> {
    Edge(EIx),
    Vertex(VIx),
}

pub struct Ubergraph<N, E, Ix> {
    // In theory I'd probably rather just do away with N
    // and make vertices merely a counter or interval tree.
    vertices: Vec<N>,
    edges: Vec<(E, Vec<EdgeMember<Ix, Ix>>)>,
}

impl<N, E> Ubergraph<N, E, usize> {
    pub fn new() -> Ubergraph<N, E, usize> {
        Ubergraph {
            vertices: Vec::new(),
            edges: Vec::new(),
        }
    }

    pub fn add_vertex(&mut self, n: N) {
        self.vertices.push(n);
    }

    pub fn add_edge(&mut self, e: E) {
        self.edges.push((e, Vec::new()));
    }

    pub fn add_node_to_edge(&mut self, idx: usize, en: EdgeMember<usize, usize>) {
        let (_, edge_nodes): &mut (E, Vec<EdgeMember<usize, usize>>) =
            &mut self.edges[idx as usize];
        edge_nodes.push(en);
    }

    pub fn edge_iter(
        &self,
    ) -> impl Iterator<Item = (&E, impl Iterator<Item = &EdgeMember<usize, usize>>)> {
        self.edges
            .iter()
            .map(|(weight, edge_vec)| (weight, edge_vec.iter()))
    }

    pub fn vert_iter(&self) -> impl Iterator<Item = &N> {
        self.vertices.iter()
    }

    fn edge_node_to_node_index(&self, edge: EdgeMember<usize, usize>) -> NodeIndex<usize> {
        match edge {
            EdgeMember::Vertex(v) => v.into(),
            EdgeMember::Edge(edge_idx) => (self.vertices.len() + edge_idx).into(),
        }
    }

    pub fn levi(&self) -> Graph<Either<&N, &E>, (), Directed, usize> {
        let mut g = Graph::<Either<&N, &E>, (), Directed, usize>::with_capacity(
            self.vertices.len() + self.edges.len(),
            0,
        );
        for vert in self.vertices.iter() {
            g.add_node(Either::Left(&vert));
        }

        for (weight, _) in self.edge_iter() {
            g.add_node(Either::Right(weight));
        }

        // must be done in a 3rd pass after adding all the vertices and edges as vertices.
        // otherwise `add_edge` could panic if we are addding an edge to an edge which
        // without an associated vertex.
        for (i, (_, nodes)) in self.edge_iter().enumerate() {
            for node in nodes {
                match *node {
                    node => {
                        g.add_edge(
                            self.edge_node_to_node_index(node),
                            self.edge_node_to_node_index(EdgeMember::Edge(i)),
                            (),
                        );
                    }
                };
            }
        }
        g
    }
}

#[derive(Copy, Clone)]
pub enum DirectedEdgeMember<VIx, EIx> {
    Member(Direction, EdgeMember<VIx, EIx>),
}

/// Lacking a formal definition, this seems appropriate given the definition of a directed hypergraph.
/// Both edges, and vertices can be either Incoming, or Outgoing.
///
/// FIXME make the types more like those of Ubergraph.
pub struct DirectedUbergraph<N, E, Ty = Directed, Ix = DirectedEdgeMember<usize, usize>> {
    vertices: Vec<N>,
    edges: Vec<(E, Vec<Ix>)>,
    ty: PhantomData<Ty>,
}

/// FIXME this hasn't really been worked on.
impl<N, E, Ty> DirectedUbergraph<N, E, Ty> {
    pub fn edge_iter(
        &self,
    ) -> impl Iterator<Item = (&E, impl Iterator<Item = &DirectedEdgeMember<usize, usize>>)> {
        self.edges
            .iter()
            .map(|(weight, edge_vec)| (weight, edge_vec.iter()))
    }
    pub fn vert_iter(&self) -> impl Iterator<Item = &N> {
        self.vertices.iter()
    }
}

#[cfg(test)]
mod tests {
    use petgraph::dot::Dot;
    use super::*;
    #[test]
    // https://arxiv.org/pdf/1704.05547.pdf
    // Example 2.
    fn example_2() {
        use EdgeMember::*;

        let mut ug = Ubergraph::new();
        let verts = [1, 2, 3];
        let edges = [
            &[Vertex(0)][..],
            &[Vertex(0), Vertex(2)][..],
            &[Vertex(0), Vertex(2), Edge(0)][..],
            &[Vertex(1), Edge(1)][..],
            &[Vertex(0), Edge(3)][..],
        ];

        for idx in verts.iter() {
            let label = format!("v{}", idx);
            ug.add_vertex(label);
        }
        for (idx, es) in edges.iter().enumerate() {
            let label = format!("e{}", idx + 1);
            ug.add_edge(label);
            for e in es.iter() {
                ug.add_node_to_edge(idx as usize, *e)
            }
        }

        let levi = ug.levi();
        println!("{:?}", Dot::with_config(&levi, &[]));
    }

    // https://arxiv.org/pdf/1704.05547.pdf
    #[test]
    fn hypergraph_example_1() {
        use EdgeMember::*;

        let mut ug = Ubergraph::new();
        let verts = [1, 2, 3, 4, 5];
        let edges = [&[Vertex(0)][..],
                     &[Vertex(0), Vertex(2)][..],
                     &[Vertex(1), Vertex(2)][..],
                     &[Vertex(0), Vertex(2), Vertex(4)][..],
        ];

        for idx in verts.iter() {
            let label = format!("v{}", idx);
            ug.add_vertex(label);
        }

        for (idx, es) in edges.iter().enumerate() {
            let label = format!("e{}", idx + 1);
            ug.add_edge(label);
            for e in es.iter() {
                ug.add_node_to_edge(idx as usize, *e)
            }
        }

        let levi = ug.levi();
        println!("{:?}", Dot::with_config(&levi, &[]));
    }
}
