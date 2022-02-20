using System.Collections;

namespace S2Geometry;

public static partial class S2ShapeUtil
{
    public class EdgeEnumerator : IEnumerator<S2Shape.Edge>, ICustomCloneable, IEquatable<EdgeEnumerator>
    {
        #region Fields, Constants

        private readonly S2ShapeIndex index_;
        private Int32 shape_id_;
        private Int32 num_edges_;
        private Int32 edge_id_;

        #endregion

        #region Constructors

        public EdgeEnumerator(S2ShapeIndex index)
        {
            index_ = index;
            shape_id_ = -1;
            num_edges_ = 0;
            edge_id_ = -1;
        }

        private EdgeEnumerator(S2ShapeIndex index, int shape_id, int num_edges, int edge_id)
        { index_ = index; shape_id_ = shape_id; num_edges_ = num_edges; edge_id_ = edge_id; }

        #endregion

        #region IEnumerator

        public S2Shape.Edge Current
        {
            get
            {
                System.Diagnostics.Debug.Assert(!Done());
                return index_.Shape(shape_id_).GetEdge(edge_id_);
            }
        }

        object IEnumerator.Current => Current;

        public bool MoveNext()
        {
            while (++edge_id_ >= num_edges_)
            {
                if (++shape_id_ >= index_.NumShapeIds()) break;
                var shape = index_.Shape(shape_id_);
                num_edges_ = (shape == null) ? 0 : shape.NumEdges();
                edge_id_ = -1;
            }
            return !Done();
        }
        public void Reset()
        {
            shape_id_ = -1;
            num_edges_ = 0;
            edge_id_ = -1;
        }
        public void Dispose() { GC.SuppressFinalize(this); }

        #endregion

        #region ICustomCloneable

        public object CustomClone() => new EdgeEnumerator(index_, shape_id_, num_edges_, edge_id_);

        #endregion

        #region EdgeEnumerator

        // Returns true if there are no more edges in the index.
        private bool Done() => shape_id_ >= index_.NumShapeIds();

        // Returns the current (shape_id, edge_id).
        public ShapeEdgeId GetShapeEdgeId() => new(shape_id_, edge_id_);

        #endregion

        #region IEquatable

        public bool Equals(EdgeEnumerator other) => index_ == other.index_ && shape_id_ == other.shape_id_ && edge_id_ == other.edge_id_;
        public override bool Equals(object obj) => obj is EdgeEnumerator ee && Equals(ee);
        public override int GetHashCode() => HashCode.Combine(index_, shape_id_, edge_id_);
        public static bool operator ==(EdgeEnumerator x, EdgeEnumerator y) => Equals(x, y);
        public static bool operator !=(EdgeEnumerator x, EdgeEnumerator y) => !Equals(x, y);

        #endregion

        #region Object

        public override string ToString() => $"(shape={shape_id_}, edge={edge_id_})";

        #endregion
    }
}
