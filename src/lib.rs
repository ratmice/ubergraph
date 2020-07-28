//! A Simple stupid attempt at building ubergraphs based on vectors.
//! While reusing a bunch of types from petgraph.
//!
//! Probably not ideal, as it means you cannot remove nodes edges, or edge members without
//! affecting the stability of indices.
//!
//! Just supports the basics of adding vertices and edges, and computing the incidence graph.

use either::Either;
use nalgebra::base::{DMatrix, Matrix, VecStorage};
use nalgebra::Dynamic;
use petgraph::graph::{Graph, NodeIndex};
use petgraph::{Directed, Direction};
use std::cmp::Ordering;
use std::iter::Iterator;
mod num_bool;

use num_bool::Bool;

/// A member of an edge
///
/// As a bit of confusing terminology sometimes I call this a node.
/// In Typical nomenclature the words `node` and `vertex` can be used interchangably.
/// When I say Vertex I mean a vertex, and when I say node, I mean a member of an edge
/// which can be either a Vertex or an Edge.
///
/// I've called this EdgeMember to avoid the confusion.
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum EdgeMember<VIx: Ord, EIx: Ord> {
    Edge(EIx),
    Vertex(VIx),
}

impl<VIx: Ord, EIx: Ord> PartialOrd for EdgeMember<VIx, EIx> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// Vertices less than Edges.
impl<VIx: Ord, EIx: Ord> Ord for EdgeMember<VIx, EIx> {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (EdgeMember::Edge(s), EdgeMember::Edge(o)) => s.cmp(o),
            (EdgeMember::Vertex(s), EdgeMember::Vertex(o)) => s.cmp(o),
            (EdgeMember::Edge(_), EdgeMember::Vertex(_)) => Ordering::Greater,
            (EdgeMember::Vertex(_), EdgeMember::Edge(_)) => Ordering::Less,
        }
    }
}

/// A recursive hypergraph structure
pub struct Ubergraph<N, E, Ix: Ord> {
    // In theory I'd probably rather just do away with N
    // and make vertices merely a counter or interval tree.
    vertices: Vec<N>,
    edges: Vec<(E, im::OrdSet<EdgeMember<Ix, Ix>>)>,
}

impl<N, E> Default for Ubergraph<N, E, usize> {
    fn default() -> Self {
        Self::new()
    }
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
        self.edges.push((e, im::OrdSet::new()));
    }

    pub fn add_node_to_edge(&mut self, idx: usize, en: EdgeMember<usize, usize>) {
        let (_, edge_nodes): &mut (E, im::OrdSet<EdgeMember<usize, usize>>) = &mut self.edges[idx];
        edge_nodes.insert(en);
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

    // FIXME
    fn matrix(
        &self,
    ) -> Matrix<
        Bool<bool>,
        nalgebra::Dynamic,
        nalgebra::Dynamic,
        VecStorage<Bool<bool>, Dynamic, Dynamic>,
    > {
        let matrix = DMatrix::<Bool<bool>>::from_element(
            self.vertices.len() + self.edges.len(),
            self.edges.len(),
            false.into(),
        );
        matrix
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
                g.add_edge(
                    self.edge_node_to_node_index(*node),
                    self.edge_node_to_node_index(EdgeMember::Edge(i)),
                    (),
                );
            }
        }
        g
    }
}

/// A member of an edge with a direction, either Incoming or Outgoing.
pub type DirectedEdgeMember<VIx, EIx> = (Direction, EdgeMember<VIx, EIx>);

impl<VIx: Ord, EIx: Ord> Into<EdgeMember<VIx, EIx>> for DirectedEdgeMember<VIx, EIx> {
    fn into(self) -> EdgeMember<VIx, EIx> {
        self.1
    }
}

/// A directed recursive ubergraph structure.
///
/// The ubergraph definition doesn't defined directed ubergraphs.
///
/// The definition given in [Directed hypergraphs and Applications](http://www.di.unipi.it/~gallo/Papers/GLNP93.ps.gz) (authors copy) [DOI](https://doi.org/10.1016/0166-218X(93)90045-P)
/// uses a enum {-1, 0, 1} matrix.
///
/// Here we use a tuple, `(petgraph::Direction, EdgeMember)`, dropping the 0
///
/// Thus a directed hyperedge would have nodes like:
/// `[(Incoming, Vertex(0)), (Outgoing, Vertex(1))]`
///
/// While a directed uberedge merely extends the Direction to edges.
/// `[(Incoming, Edge(0)) (Outgoing, Edge(1))]`
///
/// Seems like the natural thing to do, We will have to see for the Levi graph though.
pub struct DirectedUbergraph<N, E, Ix: Ord> {
    vertices: Vec<N>,
    edges: Vec<(E, Vec<DirectedEdgeMember<Ix, Ix>>)>,
}

/// FIXME this hasn't really been worked on.
impl<N, E> DirectedUbergraph<N, E, usize> {
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
    use super::*;
    use insta;
    use petgraph::dot::Dot;

    fn test_helper<T: std::fmt::Display>(
        vertices: &[T],
        edges: &[&[EdgeMember<usize, usize>]],
    ) -> Ubergraph<String, String, usize> {
        let mut ug = Ubergraph::new();
        for idx in vertices.iter() {
            let label = format!("v{}", idx);
            ug.add_vertex(label);
        }
        for (idx, es) in edges.iter().enumerate() {
            let label = format!("e{}", idx + 1);
            ug.add_edge(label);
            for e in es.iter() {
                ug.add_node_to_edge(idx, *e)
            }
        }
        ug
    }

    // https://arxiv.org/pdf/1704.05547.pdf
    // Example 2.
    #[test]
    fn example_2() {
        use EdgeMember::*;
        let verts = [1, 2, 3];
        let edges = [
            &[Vertex(0)][..],
            &[Vertex(0), Vertex(2)][..],
            &[Vertex(0), Vertex(2), Edge(0)][..],
            &[Vertex(1), Edge(1)][..],
            &[Vertex(0), Edge(3)][..],
        ];

        let ug = test_helper(&verts, &edges);
        // For some reason this output looks nicer than assert_debug_snapshot!().
        let dot_str = format!("{:?}", Dot::with_config(&ug.levi(), Default::default()));
        insta::assert_snapshot!(dot_str);
    }

    // https://arxiv.org/pdf/1704.05547.pdf
    #[test]
    fn hypergraph_example_1() {
        use EdgeMember::*;
        let verts = [1, 2, 3, 4, 5];
        let edges = [
            &[Vertex(0)][..],
            &[Vertex(0), Vertex(2)][..],
            &[Vertex(1), Vertex(2)][..],
            &[Vertex(0), Vertex(2), Vertex(4)][..],
        ];
        let ug = test_helper(&verts, &edges);
        // For some reason this output looks nicer than assert_debug_snapshot!().
        let dot_str = format!("{:?}", Dot::with_config(&ug.levi(), Default::default()));
        insta::assert_snapshot!(dot_str);
    }
}
