using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;

namespace S2Geometry
{
    /// <summary>
    /// Given a sequence of S2Points assumed to be the center of level-k cells,
    /// compresses it into a stream using the following method:
    /// - decompose the points into (face, si, ti) tuples (see s2coords.h)
    /// - run-length encode the faces, combining face number and count into a
    ///     varInt32.  See the Faces class in s2point_compression.cc.
    /// - right shift the (si, ti) to remove the part that's ant for all cells
    ///     of level-k.  The result is called the (pi, qi) space.
    /// - 2nd derivative encode the pi and qi sequences (linear prediction)
    /// - zig-zag encode all derivative values but the first, which cannot be
    ///     negative
    /// - interleave the zig-zag encoded values
    /// - encode the first interleaved value in a fixed length encoding
    ///     (varint would make this value larger)
    /// - encode the remaining interleaved values as varint64s, as the
    ///     derivative encoding should make the values small.
    /// In addition, provides a lossless method to compress a sequence of points even
    /// if some points are not the center of level-k cells. These points are stored
    /// exactly, using 3 double precision values, after the above encoded string,
    /// together with their index in the sequence (this leads to some redundancy - it
    /// is expected that only a small fraction of the points are not cell centers).
    /// 
    /// Require that the encoder was constructed with the no-arg constructor, as
    /// Ensure() will be called to allocate space.
    /// 
    /// 
    /// To encode leaf cells, this requires 8 bytes for the first vertex plus
    /// an average of 3.8 bytes for each additional vertex, when computed on
    /// Google's geographic repository.
    /// </summary>
    public static class S2PointCompression
    {
        private const int kDerivativeEncodingOrder = 2;

        /// <summary>
        /// The XYZ and face,si,ti coordinates of an S2Point and, if this point is equal
        /// to the center of an S2Cell, the level of this cell (-1 otherwise).
        /// </summary>
        public record S2XYZFaceSiTi(S2Point XYZ, int Face, uint Si, uint Ti, int CellLevel);

        /// <summary>
        /// Pair of face number and count for run-length encoding.
        /// </summary>
        private class FaceRun
        {
            public readonly int Face;
            public int Count;

            public FaceRun()
            {
                Face = -1;
                Count = 0;
            }
            public FaceRun(int initial_face, int initial_count)
            {
                Face = initial_face;
                Count = initial_count;
            }

            /// <summary>
            /// Encodes each face as a varint64 with value kNumFaces * count + face.
            /// 21 faces can fit in a single byte.  Varint64 is used so that 4G faces
            /// can be encoded instead of just 4G / 6 = ~700M.
            /// </summary>
            public void Encode(Encoder encoder)
            {
                encoder.Ensure(Encoder.kVarintMax64);

                // It isn't necessary to encode the number of faces left for the last run,
                // but since this would only help if there were more than 21 faces, it will
                // be a small overall savings, much smaller than the bound encoding.
                encoder.PutVarUInt64((UInt64)(S2CellId.kNumFaces * (Int64)Count + Face));
                Assert.True(encoder.Avail() >= 0);
            }

            public static bool Decode(Decoder decoder, out FaceRun fr)
            {
                fr = null;
                if (!decoder.TryGetVarUInt64(out var face_and_count)) return false;

                var face = (int)(face_and_count % S2CellId.kNumFaces);
                // Make sure large counts don't wrap on malicious or random input.
                UInt64 count64 = face_and_count / S2CellId.kNumFaces;
                var count = (int)count64;

                var res = count > 0 && (UInt64)count == count64;
                if (!res) return false;

                fr = new FaceRun(face, count);
                return true;
            }
        }

        /// <summary>
        /// Run-length encoder/decoder for face numbers.
        /// </summary>
        private class Faces
        {
            public class Enumerator : IEnumerator<int>
            {
                /// <summary>
                /// The faces_ vector of the Faces object for which this is an iterator.
                /// </summary>
                private readonly List<FaceRun> faces_;

                /// <summary>
                /// The index that the next face will come from.
                /// </summary>
                private int face_index_;

                // Number of faces already consumed for face_index_.
                private int num_faces_used_for_index_;

                public Enumerator(Faces faces)
                {
                    faces_ = faces.faces_;
                    face_index_ = 0;
                    num_faces_used_for_index_ = 0;
                }

                public int Current => faces_[face_index_].Face;

                object IEnumerator.Current => Current;

                public void Dispose() => GC.SuppressFinalize(this);

                public bool MoveNext()
                {
                    Assert.True(faces_.Count != face_index_);
                    Assert.True(num_faces_used_for_index_ <= faces_[face_index_].Count);
                    if (num_faces_used_for_index_ == faces_[face_index_].Count)
                    {
                        face_index_++;
                        num_faces_used_for_index_ = 0;
                    }

                    num_faces_used_for_index_++;
                    return face_index_ < faces_.Count;
                }

                public void Reset() => face_index_ = num_faces_used_for_index_ = 0;
            }

            /// <summary>
            /// Add the face to the list of face runs, combining with the last if
            /// possible.
            /// </summary>
            public void AddFace(int face)
            {
                if (faces_.Any() && faces_.Last().Face == face)
                {
                    faces_.Last().Count++;
                }
                else
                {
                    faces_.Add(new FaceRun(face, 1));
                }
            }

            // Encodes the faces to encoder.
            public void Encode(Encoder encoder)
            {
                foreach (var face_run in faces_)
                    face_run.Encode(encoder);
            }

            // Decodes the faces, returning true on success.
            public bool Decode(int num_vertices, Decoder decoder)
            {
                int num_faces_parsed = 0;
                while (num_faces_parsed < num_vertices)
                {
                    if (!FaceRun.Decode(decoder, out var face_run)) return false;

                    faces_.Add(face_run);
                    num_faces_parsed += face_run.Count;
                }

                return true;
            }

            public Enumerator GetEnumerator()
            {
                return new Enumerator(this);
            }

            /// <summary>
            /// Run-length encoded list of faces.
            /// </summary>
            private readonly List<FaceRun> faces_= new();
        }

        /// <summary>
        /// Unused function (for documentation purposes only).
        /// </summary>
#pragma warning disable IDE0051 // Quitar miembros privados no utilizados
        private static int STtoPiQi(double s, int level)
#pragma warning restore IDE0051 // Quitar miembros privados no utilizados
        {
            // We introduce a new coordinate system (pi, qi), which is (si, ti)
            // with the bits that are constant for cells of that level shifted
            // off to the right.
            // si = round(s * 2^31)
            // pi = si >> (31 - level)
            //    = floor(s * 2^level)
            // If the point has been snapped to the level, the bits that are
            // shifted off will be a 1 in the msb, then 0s after that, so the
            // fractional part discarded by the cast is (close to) 0.5.
            return (int)(s * (1 << level));
        }

        private static int SiTitoPiQi(uint si, int level)
        {
            // See STtoPiQi for the definition of the PiQi coordinate system.
            //
            // EncodeFirstPointFixedLength encodes the return value using "level" bits,
            // so we clamp "si" to the range [0, 2**level - 1] before trying to encode
            // it.  This is okay because if si == kMaxSiTi, then it is not a cell center
            // anyway and will be encoded separately as an "off-center" point.
            si = Math.Min(si, S2Constants.kMaxSiTi - 1);
            return (int)(si >> (S2Constants.kMaxCellLevel + 1 - level));
        }

        private static double PiQitoST(int pi, int level)
        {
            // We want to recover the position at the center of the cell.  If the point
            // was snapped to the center of the cell, then modf(s * 2^level) == 0.5.
            // Inverting STtoPiQi gives:
            // s = (pi + 0.5) / 2^level.
            return (pi + 0.5) / (1 << level);
        }

        private static S2Point FacePiQitoXYZ(int face, int pi, int qi, int level)
        {
            return S2Coords.FaceUVtoXYZ(face,
                S2Coords.STtoUV(PiQitoST(pi, level)),
                S2Coords.STtoUV(PiQitoST(qi, level))).Normalized;
        }

        private static void EncodeFirstPointFixedLength((int, int) vertex_pi_qi, int level, NthDerivativeCoder pi_coder, NthDerivativeCoder qi_coder, Encoder encoder)
        {
            // Do not ZigZagEncode the first point, since it cannot be negative.
            UInt32 pi = (uint)pi_coder.Encode(vertex_pi_qi.Item1);
            UInt32 qi = (uint)qi_coder.Encode(vertex_pi_qi.Item2);
            // Interleave to reduce overhead from two partial bytes to one.
            UInt64 interleaved_pi_qi = BitsInterleave.InterleaveUInt32(pi, qi);

            // Convert to little endian for architecture independence.
            var little_endian_interleaved_pi_qi = BitConverter.GetBytes(interleaved_pi_qi);

            int bytes_required = (level + 7) / 8 * 2;
            Assert.True(bytes_required <= 8);
            encoder.Ensure(bytes_required);
            encoder.PutN(little_endian_interleaved_pi_qi, bytes_required);
            Assert.True(encoder.Avail() >= 0);
        }

#pragma warning disable IDE0060 // Quitar el parámetro no utilizado
        private static void EncodePointCompressed((int, int) vertex_pi_qi, int level, NthDerivativeCoder pi_coder, NthDerivativeCoder qi_coder, Encoder encoder)
#pragma warning restore IDE0060 // Quitar el parámetro no utilizado
        {
            // ZigZagEncode, as varint requires the maximum number of bytes for
            // negative numbers.
            UInt32 zig_zag_encoded_deriv_pi = Transforms.ZigZagEncode(pi_coder.Encode(vertex_pi_qi.Item1));
            UInt32 zig_zag_encoded_deriv_qi = Transforms.ZigZagEncode(qi_coder.Encode(vertex_pi_qi.Item2));
            // Interleave to reduce overhead from two partial bytes to one.
            UInt64 interleaved_zig_zag_encoded_derivs = BitsInterleave.InterleaveUInt32(zig_zag_encoded_deriv_pi, zig_zag_encoded_deriv_qi);

            encoder.Ensure(Encoder.kVarintMax64);
            encoder.PutVarUInt64(interleaved_zig_zag_encoded_derivs); // TODO: check if it would be better to use signed int and remove zigzag encoding here, modify "DecodePointCompressed" accordingly
            Assert.True(encoder.Avail() >= 0);
        }

        private static void EncodePointsCompressed(List<(int, int)> vertices_pi_qi, int level, Encoder encoder)
        {
            var pi_coder = new NthDerivativeCoder(kDerivativeEncodingOrder);
            var qi_coder = new NthDerivativeCoder(kDerivativeEncodingOrder);
            for (int i = 0; i < vertices_pi_qi.Count; ++i)
            {
                if (i == 0)
                {
                    // The first point will be just the (pi, qi) coordinates
                    // of the S2Point.  NthDerivativeCoder will not save anything
                    // in that case, so we encode in fixed format rather than varint
                    // to avoid the varint overhead.
                    EncodeFirstPointFixedLength(vertices_pi_qi[i], level, pi_coder, qi_coder, encoder);
                }
                else
                {
                    EncodePointCompressed(vertices_pi_qi[i], level, pi_coder, qi_coder, encoder);
                }
            }

            Assert.True(encoder.Avail() >= 0);
        }

        private static bool DecodeFirstPointFixedLength(Decoder decoder, int level, NthDerivativeCoder pi_coder, NthDerivativeCoder qi_coder, out (int, int) vertex_pi_qi)
        {
            vertex_pi_qi = default;
            int bytes_required = (level + 7) / 8 * 2;
            if (decoder.Avail() < bytes_required) return false;
            var little_endian_interleaved_pi_qi = new byte[8];

            decoder.GetN(little_endian_interleaved_pi_qi, 0, bytes_required);

            UInt64 interleaved_pi_qi = BitConverter.ToUInt64(little_endian_interleaved_pi_qi, 0);

            BitsInterleave.DeinterleaveUInt32(interleaved_pi_qi, out var pi, out var qi);

            var t1 = pi_coder.Decode((int)pi);
            var t2 = qi_coder.Decode((int)qi);
            vertex_pi_qi = (t1, t2);
            return true;
        }

#pragma warning disable IDE0060 // Quitar el parámetro no utilizado
        public static bool DecodePointCompressed(Decoder decoder, int level, NthDerivativeCoder pi_coder, NthDerivativeCoder qi_coder, out (int, int) vertex_pi_qi)
#pragma warning restore IDE0060 // Quitar el parámetro no utilizado
        {
            if (!decoder.TryGetVarUInt64(out var interleaved_zig_zag_encoded_deriv_pi_qi))
            {
                vertex_pi_qi = default;
                return false;
            }

            BitsInterleave.DeinterleaveUInt32(interleaved_zig_zag_encoded_deriv_pi_qi,
                                          out var zig_zag_encoded_deriv_pi,
                                          out var zig_zag_encoded_deriv_qi);

            var t1 = pi_coder.Decode(Transforms.ZigZagDecode(zig_zag_encoded_deriv_pi));
            var t2 = qi_coder.Decode(Transforms.ZigZagDecode(zig_zag_encoded_deriv_qi));
            vertex_pi_qi = (t1, t2);
            return true;
        }

        /// <summary>
        /// Encode the points in the encoder, using an optimized compressed format for
        /// points at the center of a cell at 'level', plus 3 double values for the
        /// others.
        /// </summary>
        public static void S2EncodePointsCompressed(S2XYZFaceSiTi[] points, int level, Encoder encoder)
        {
            var vertices_pi_qi = new List<(int, int)>(points.Length);
            var off_center = new List<int>();
            var faces = new Faces();
            for (int i = 0; i < points.Length; ++i)
            {
                faces.AddFace(points[i].Face);
                vertices_pi_qi.Add((
                    SiTitoPiQi(points[i].Si, level),
                    SiTitoPiQi(points[i].Ti, level)));
                if (points[i].CellLevel != level)
                {
                    off_center.Add(i);
                }
            }
            faces.Encode(encoder);
            EncodePointsCompressed(vertices_pi_qi, level, encoder);
            var num_off_center = off_center.Count;
            encoder.Ensure(Encoder.kVarintMax32 + (Encoder.kVarintMax32 + Marshal.SizeOf(typeof(S2Point))) * num_off_center);
            encoder.PutVarUInt32((uint)num_off_center);
            Assert.True(encoder.Avail() >= 0);
            foreach (var index in off_center)
            {
                encoder.PutVarUInt32((uint)index);
                encoder.PutPoints(new []{ points[index].XYZ });
                Assert.True(encoder.Avail() >= 0);
            }
        }

        /// <summary>
        /// Decode points encoded with S2EncodePointsCompressed. Requires that the
        /// level is the level that was used in S2EncodePointsCompressed. Ensures
        /// that the decoded points equal the encoded points. Returns true on success.
        /// </summary>
        public static bool S2DecodePointsCompressed(Decoder decoder, int level, S2Point[] points, int offset)
        {
            var faces = new Faces();
            var len = points.Length - offset;
            if (!faces.Decode(len, decoder))
            {
                return false;
            }

            var pi_coder = new NthDerivativeCoder(kDerivativeEncodingOrder);
            var qi_coder = new NthDerivativeCoder(kDerivativeEncodingOrder);
            var faces_iterator = faces.GetEnumerator();
            for (int i = 0; i < len; ++i)
            {
                (int, int) vertex_pi_qi;
                if (i == 0)
                {
                    if (!DecodeFirstPointFixedLength(decoder, level, pi_coder, qi_coder, out vertex_pi_qi))
                    {
                        return false;
                    }
                }
                else
                {
                    if (!DecodePointCompressed(decoder, level, pi_coder, qi_coder, out vertex_pi_qi))
                    {
                        return false;
                    }
                }
                
                faces_iterator.MoveNext();
                int face = faces_iterator.Current;
                points[offset + i] = FacePiQitoXYZ(face, vertex_pi_qi.Item1, vertex_pi_qi.Item2, level);
            }

            if (!decoder.TryGetVarUInt32(out var num_off_center) || num_off_center > len)
            {
                return false;
            }
            for (int i = 0; i < num_off_center; ++i)
            {
                if (!decoder.TryGetVarUInt32(out var index) || index >= len)
                {
                    return false;
                }
                if (decoder.Avail() < Marshal.SizeOf(typeof(S2Point))) return false;
                decoder.GetPoints(points, offset + (int)index, 1);
            }
            return true;
        }
    }
}
