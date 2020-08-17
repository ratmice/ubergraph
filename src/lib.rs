//! A Simple stupid attempt at building ubergraphs based on vectors.
//! While reusing a bunch of types from petgraph.
//!
//! Probably not ideal, as it means you cannot remove nodes edges, or edge members without
//! affecting the stability of indices.
//!
//! Just supports the basics of adding vertices and edges, and computing the incidence graph.

use either::Either;
use itertools::Itertools;
use nalgebra::base::{Matrix, MatrixMN, VecStorage};
use nalgebra::Dynamic;
use petgraph::graph::{Graph, NodeIndex};
use petgraph::{Directed, Direction};
use std::cmp::Ordering;
use std::iter::Iterator;

mod indices;
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
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum EdgeMember<VIx: Ord, EIx: Ord> {
    Vertex(VIx),
    Edge(EIx),
}

impl<VIx: Ord, EIx: Ord> EdgeMember<VIx, EIx> {
    pub fn is_edge(self) -> bool {
        if let EdgeMember::Edge(_) = self {
            true
        } else {
            false
        }
    }
    pub fn is_vertex(self) -> bool {
        if let EdgeMember::Vertex(_) = self {
            true
        } else {
            false
        }
    }
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
/// Currently not well-founded, in that the construction allows you to introduce cycles.
pub struct Ubergraph<V, W, Ix: Ord> {
    vertices: Vec<V>,
    edges: Vec<(W, im::OrdSet<EdgeMember<Ix, Ix>>)>,
}

#[derive(Copy, Clone, PartialEq, Eq, Debug, Ord, PartialOrd)]
pub enum HyperEdgeMember<VIx: Ord> {
    Vertex(VIx),
}

///
/// This is not ideal, as the ubergraph.edges can't be dropped,
/// but we don't really want to clone the vertices vector
/// the actual ownership and field types here are pretty certainly
/// assured to change.
pub struct HypergraphRepresentation<'a, V, W, Ix: Ord> {
    ubergraph: &'a Ubergraph<V, W, Ix>,
    edges: Vec<im::OrdSet<HyperEdgeMember<Ix>>>,
}

impl<V, W> HypergraphRepresentation<'_, V, W, usize> {
    pub fn matrix(
        &self,
    ) -> Matrix<
        Bool<bool>,
        nalgebra::Dynamic,
        nalgebra::Dynamic,
        VecStorage<Bool<bool>, Dynamic, Dynamic>,
    > {
        MatrixMN::<Bool<bool>, nalgebra::Dynamic, nalgebra::Dynamic>::from_iterator(
            self.ubergraph.vertices.len(),
            self.edges.len(),
            self.edges.iter().flat_map(|edge_set| {
                let mut pos = 0;
                edge_set
                    .iter()
                    .flat_map(move |edge_mem| {
                        let idx = match *edge_mem {
                            HyperEdgeMember::Vertex(x) => x,
                        };
                        let it = std::iter::repeat(false.into())
                            .take(idx - pos)
                            .chain(std::iter::once(true.into()));
                        pos = idx + 1;
                        it
                    })
                    .pad_using(self.ubergraph.vertices.len(), |_| false.into())
            }),
        )
    }
    pub fn levi(&self) -> Graph<Either<&V, &W>, (), Directed, usize> {
        let mut g = Graph::<Either<&V, &W>, (), Directed, usize>::with_capacity(
            self.ubergraph.vertices.len() + self.edges.len(),
            0,
        );

        for vert in self.ubergraph.vertices.iter() {
            g.add_node(Either::Left(&vert));
        }

        // unlike for ubergraph this is fine with 2 passes,
        // since we've already added all the vertices, there no possibility to run across an edge
        // directed at anything which hasn't been added.
        for (eidx, edge_set) in self.edges.iter().enumerate() {
            g.add_node(Either::Right(&self.ubergraph.edges[eidx].0));
            for edge_mem in edge_set.iter() {
                match edge_mem {
                    HyperEdgeMember::Vertex(vidx) => g.add_edge(
                        (*vidx).into(),
                        (self.ubergraph.vertices.len() + eidx).into(),
                        (),
                    ),
                };
            }
        }
        g
    }
}

impl<V, W> Default for Ubergraph<V, W, usize> {
    fn default() -> Self {
        Self::new()
    }
}

impl<V, W> Ubergraph<V, W, usize> {
    pub fn new() -> Ubergraph<V, W, usize> {
        Ubergraph {
            vertices: Vec::new(),
            edges: Vec::new(),
        }
    }

    pub fn add_vertex(&mut self, n: V) {
        self.vertices.push(n);
    }

    pub fn add_edge(&mut self, w: W) {
        self.edges.push((w, im::OrdSet::new()));
    }

    pub fn add_node_to_edge(&mut self, idx: usize, en: EdgeMember<usize, usize>) {
        let (_, edge_nodes): &mut (W, im::OrdSet<EdgeMember<usize, usize>>) = &mut self.edges[idx];
        edge_nodes.insert(en);
    }

    pub fn edge_iter(
        &self,
    ) -> impl Iterator<Item = (&W, impl Iterator<Item = &EdgeMember<usize, usize>>)> {
        self.edges
            .iter()
            .map(|(weight, edge_vec)| (weight, edge_vec.iter()))
    }

    pub fn vert_iter(&self) -> impl Iterator<Item = &V> {
        self.vertices.iter()
    }

    fn edge_node_to_node_index(&self, edge: EdgeMember<usize, usize>) -> NodeIndex<usize> {
        self.column_index(edge).into()
    }

    fn column_index(&self, edge: EdgeMember<usize, usize>) -> usize {
        match edge {
            EdgeMember::Vertex(v) => v,
            EdgeMember::Edge(edge_idx) => (self.vertices.len() + edge_idx),
        }
    }

    /// nalgebra probably *isn't* the right matrix to be using for a boolean matrix,
    /// It seems a bit overkill, and it requires us to convert bool into a numeric type.
    /// However, nalgebra does work with a non-square matrix,
    /// Which is why i'm using it anyway for now.
    /// The returned matrix will be N*M+M in size
    pub fn matrix(
        &self,
    ) -> Matrix<
        Bool<bool>,
        nalgebra::Dynamic,
        nalgebra::Dynamic,
        VecStorage<Bool<bool>, Dynamic, Dynamic>,
    > {
        MatrixMN::<Bool<bool>, nalgebra::Dynamic, nalgebra::Dynamic>::from_iterator(
            self.vertices.len() + self.edges.len(),
            self.edges.len(),
            self.edges.iter().flat_map(|(_, edge_set)| {
                let mut pos = 0;
                edge_set
                    .iter()
                    .flat_map(move |edge_mem| {
                        let idx = self.column_index(*edge_mem);
                        let it = std::iter::repeat(false.into())
                            .take(idx - pos)
                            .chain(std::iter::once(true.into()));
                        pos = idx + 1;
                        it
                    })
                    .pad_using(self.vertices.len() + self.edges.len(), |_| false.into())
            }),
        )
    }

    pub fn levi(&self) -> Graph<Either<&V, &W>, (), Directed, usize> {
        let mut g = Graph::<Either<&V, &W>, (), Directed, usize>::with_capacity(
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
        for (eidx, (_, edge_set)) in self.edge_iter().enumerate() {
            for edge_mem in edge_set {
                g.add_edge(
                    self.edge_node_to_node_index(*edge_mem),
                    self.edge_node_to_node_index(EdgeMember::Edge(eidx)),
                    (),
                );
            }
        }
        g
    }

    fn hypergraph_expansion(
        &self,
        edge_set: &im::OrdSet<EdgeMember<usize, usize>>,
    ) -> im::OrdSet<EdgeMember<usize, usize>> {
        edge_set
            .iter()
            .flat_map(|edge_mem| match *edge_mem {
                EdgeMember::Edge(x) => self.hypergraph_expansion(&self.edges[x].1),
                vert => im::OrdSet::unit(vert),
            })
            .collect()
    }

    /// Returns an ubergraph, with no EdgeMembers that are EdgeMember::Edge,
    /// Edge to Edge links are collapsed recursively into Vertices.
    /// In this sense it is still an ubergraph with N + M * M matrix, but the + M rows are all
    /// zero.
    pub fn hypergraph(&self) -> HypergraphRepresentation<'_, V, W, usize> {
        HypergraphRepresentation {
            ubergraph: self,
            edges: self
                .edges
                .iter()
                .map(|(_, edge_set)| {
                    self.hypergraph_expansion(edge_set)
                        .iter()
                        .map(|edge_mem| match edge_mem {
                            EdgeMember::Vertex(v) => HyperEdgeMember::Vertex(*v),
                            _ => unreachable!(),
                        })
                        .collect()
                })
                .collect(),
        }
    }

    fn relative_index(&self, edge: EdgeMember<usize, usize>) -> usize {
        match edge {
            EdgeMember::Vertex(v) => v,
            EdgeMember::Edge(edge_idx) => edge_idx,
        }
    }

    // Probably need to hand roll this with loop labels
    // instead of recursion to avoid overflowing the stack.
    fn edge_depth(&self, edge: &EdgeMember<usize, usize>, k: usize) -> usize {
        assert!(EdgeMember::is_edge(*edge));
        let idx = self.relative_index(*edge);
        let mut max_depth = k;
        for edge in self.edges[idx].1.iter().filter(|edge| edge.is_edge()) {
            let depth = self.edge_depth(edge, k + 1);
            if depth > max_depth {
                max_depth = depth
            }
        }
        max_depth
    }

    pub fn depth(&self) -> usize {
        let mut k = 0;
        assert!(self.edges.len() > 0);
        // FIXME We shouldn't traverse edges which were already visited.
        for (idx, _) in self.edges.iter().enumerate() {
            let depth = self.edge_depth(&EdgeMember::Edge(idx), 0);
            if depth > k {
                k = depth
            }
        }
        k
    }
}

/// A member of an edge with a direction, either Incoming or Outgoing.
///
/// The following seems like the natural thing to do, however can lead
/// not well-founded graphs, consider the graph e1 = {(Incoming, e2), (Outgoing, e2)}
pub type DirectedEdgeMember<VIx, EIx> = (Direction, EdgeMember<VIx, EIx>);

impl<VIx: Ord, EIx: Ord> Into<EdgeMember<VIx, EIx>> for DirectedEdgeMember<VIx, EIx> {
    fn into(self) -> EdgeMember<VIx, EIx> {
        self.1
    }
}

/// Maybe if we throw more specializations in here a nice generalization will appear.
/// Eventually we'll need to clean all this up.
///
/// A second definition of edge member for edges with directed vertices,
/// either `Incoming` or `Outgoing`.
///
/// Where the levy graph may only contains Edge -> Vertex cycles.
/// As such all paths through edges, not traversing vertices remain acyclic.
///
/// This seems like a better definition which should be sufficient for my purposes...
pub enum VertexDirectedEdgeMember<VIx: Ord, EIx: Ord> {
    Vertex(Direction, VIx),
    Edge(EIx),
}

/// A directed recursive ubergraph structure.
///
/// The ubergraph definition doesn't defined directed ubergraphs.
///
/// The definition somewhat related to [Directed hypergraphs and Applications](http://www.di.unipi.it/~gallo/Papers/GLNP93.ps.gz) (authors copy) [DOI](https://doi.org/10.1016/0166-218X(93)90045-P)
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

/// Arrows point from source to target.
/// Basically The definition from [Directed hypergraphs and Applications](http://www.di.unipi.it/~gallo/Papers/GLNP93.ps.gz) (authors copy) [DOI](https://doi.org/10.1016/0166-218X(93)90045-P)
pub struct DirectedHyperEdge<Ix> {
    source: im::OrdSet<Ix>,
    target: im::OrdSet<Ix>,
}

pub struct GalloDirectedUbergraph<N, E, Ix: Ord> {
    vertices: Vec<N>,
    edges: (E, DirectedHyperEdge<Ix>),
}

// In theory the above definition covers these below, but
// these are more specialized.
pub struct FArc<Ix> {
    source: Ix,
    target: im::OrdSet<Ix>,
}

pub struct BArc<Ix> {
    source: im::OrdSet<Ix>,
    target: Ix,
}

enum BFArc<Ix> {
    FArc(FArc<Ix>),
    BArc(BArc<Ix>),
}

pub struct BGraph<N, E, Ix: Ord> {
    vertices: Vec<N>,
    edges: (E, BArc<Ix>),
}

pub struct FGraph<N, E, Ix: Ord> {
    vertices: Vec<N>,
    edges: (E, FArc<Ix>),
}

pub struct BFGraph<N, E, Ix: Ord> {
    vertices: Vec<N>,
    edges: (E, BFArc<Ix>),
}

#[cfg(test)]
mod tests {
    use super::*;
    use insta;
    use petgraph::dot::Dot;
    use EdgeMember::*;

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

    struct TestGraph<'a, T> {
        verts: &'a [T],
        edges: &'a [&'a [EdgeMember<usize, usize>]],
    }

    // Example 2 from https://arxiv.org/pdf/1704.05547.pdf
    const EXAMPLE2: TestGraph<u32> = TestGraph {
        verts: &[1, 2, 3],
        edges: &[
            &[Vertex(0)],
            &[Vertex(0), Vertex(2)],
            &[Vertex(0), Vertex(2), Edge(0)],
            &[Vertex(1), Edge(1)],
            &[Vertex(0), Edge(3)],
        ],
    };

    #[test]
    fn ubergraph_example_2() {
        let ug = test_helper(&EXAMPLE2.verts, &EXAMPLE2.edges);
        // For some reason this output looks nicer than assert_debug_snapshot!().
        let dot_str = format!("{:?}", Dot::with_config(&ug.levi(), Default::default()));
        insta::assert_snapshot!(dot_str);
    }

    #[test]
    #[rustfmt::skip]
    fn example2_matrix() {
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);
        let matrix = ug.matrix();

        assert_eq!(
             matrix.as_slice(),
            vec![ true, false, false, false, false, false, false, false,
                  true, false, true, false, false, false, false, false,
                  true, false, true, true, false, false, false, false,
                  false, true, false, false, true, false, false, false,
                  true, false, false, false, false, false, true, false ]
            .iter()
            .map(|b| Bool::from(*b))
            .collect::<Vec<Bool<bool>>>()
            .as_slice()
        )
    }

    #[test]
    #[rustfmt::skip]
    fn example2_to_hypergraph_matrix() {
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);
        assert_eq!(ug.hypergraph().matrix().as_slice(),
            vec![ true, false, false,
                  true, false, true,
                  true, false, true,
                  true, true,  true,
                  true, true,  true,  ]
            .iter()
            .map(|b| Bool::from(*b))
            .collect::<Vec<Bool<bool>>>()
            .as_slice()
        );
    }

    #[test]
    fn ubergraph_depth() {
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);

        assert_eq!(ug.depth(), 2);
    }
    #[test]
    fn ubergraph_edge_depth() {
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);

        assert_eq!(ug.edge_depth(&EdgeMember::Edge(1), 0), 0);
        assert_eq!(ug.edge_depth(&EdgeMember::Edge(2), 0), 1);
        assert_eq!(ug.edge_depth(&EdgeMember::Edge(4), 0), 2);
    }

    #[test]
    fn ubergraph_to_hypergraph_to_levi_graph() {
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);
        let hg = ug.hypergraph();
        let levi = hg.levi();

        let dot_str = format!("{:?}", Dot::with_config(&levi, Default::default()));
        insta::assert_snapshot!(dot_str);
    }

    #[test]
    fn check_matrix_impl() {
        // The matrix constructor *should* be more efficient, but this one is more obviously
        // correct.
        // TODO check this with critereon...
        let ug = test_helper(EXAMPLE2.verts, EXAMPLE2.edges);
        assert_eq!(
            ug.matrix().data.as_vec(),
            &ug.edges
                .iter()
                .flat_map(|(_, edge_set)| {
                    (0..ug.vertices.len())
                        .map(move |idx| edge_set.contains(&EdgeMember::Vertex(idx)).into())
                        .chain(
                            (0..ug.edges.len())
                                .map(move |idx| edge_set.contains(&EdgeMember::Edge(idx)).into()),
                        )
                })
                .collect::<Vec<Bool<bool>>>()
        )
    }

    // Example 1 from https://arxiv.org/pdf/1704.05547.pdf
    const EXAMPLE1: TestGraph<u32> = TestGraph {
        verts: &[1, 2, 3, 4, 5],
        edges: &[
            &[Vertex(0)],
            &[Vertex(0), Vertex(2)],
            &[Vertex(1), Vertex(2)],
            &[Vertex(0), Vertex(2), Vertex(4)],
        ],
    };

    #[test]
    fn hypergraph_example_1() {
        let ug = test_helper(&EXAMPLE1.verts, &EXAMPLE1.edges);
        // For some reason this output looks nicer than assert_debug_snapshot!().
        let dot_str = format!("{:?}", Dot::with_config(&ug.levi(), Default::default()));
        insta::assert_snapshot!(dot_str);
    }

    #[test]
    fn hypergraph_depth() {
        let ug = test_helper(&EXAMPLE1.verts, &EXAMPLE1.edges);
        // All hypergraphs should be depth 0.
        assert_eq!(ug.depth(), 0);
    }
}
