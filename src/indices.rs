struct VertexIndex<VT>(VT);
struct EdgeIndex<ET>(ET);
struct LastVertexIndex<VT>(VT);
struct MatrixIndex<T>(T);

use std::ops::Add;

impl<ET,VT> std::ops::Add<EdgeIndex<ET>> for LastVertexIndex<VT>
  where VT: std::ops::Add<ET>,
{
  type Output = MatrixIndex<<VT as std::ops::Add<ET>>::Output>;

  fn add(self: LastVertexIndex<VT>, other: EdgeIndex<ET>) -> MatrixIndex<<VT as Add<ET>>::Output> {
        MatrixIndex(self.0 + other.0)
  }
}

