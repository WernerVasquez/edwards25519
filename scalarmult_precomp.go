package edwards25519

// precompTable is used to quick scalar multiplication when the same Point is reused
type precompTable [32]affineLookupTable

// pointTablePrecomp maps the canonical 32-byte encoding of a Point to its precomputed table
// this is used to store a previously precomputed Point so that it does not need to be recalculated
var pointTablePrecomp = make(map[[32]byte]precompTable)

// PrecompPoint returns *precompTable
// It is precomputed the first time it's used for each Point p.
// It is the only way to create the *precompTable
// calling this has appreciable overhead even if it is just reusing previous computed results
func PrecompPoint(p *Point) *precompTable {

	pBytes := (*[32]byte)(p.Bytes())
	var v precompTable
	found := false

	if v, found = pointTablePrecomp[*pBytes]; !found {
		pCopy := new(Point).Set(p)
		for i := 0; i < 32; i++ {
			v[i].FromP3(pCopy)
			for j := 0; j < 8; j++ {
				pCopy.Add(pCopy, pCopy)
			}
		}
		pointTablePrecomp[*pBytes] = v
	}

	return &v
}

// ScalarMult returns x * P.
//
// The scalar multiplication is done in constant time.
// This is as fast as ScalarBaseMult
func (pointTable *precompTable) ScalarMult(x *Scalar) *Point {
	//pointTable := pointTable(p)

	// Write x = sum(x_i * 16^i) so  x*B = sum( B*x_i*16^i )
	// as described in the Ed25519 paper
	//
	// Group even and odd coefficients
	// x*B     = x_0*16^0*B + x_2*16^2*B + ... + x_62*16^62*B
	//         + x_1*16^1*B + x_3*16^3*B + ... + x_63*16^63*B
	// x*B     = x_0*16^0*B + x_2*16^2*B + ... + x_62*16^62*B
	//    + 16*( x_1*16^0*B + x_3*16^2*B + ... + x_63*16^62*B)
	//
	// We use a lookup table for each i to get x_i*16^(2*i)*B
	// and do four doublings to multiply by 16.
	digits := x.signedRadix16()

	multiple := &affineCached{}
	tmp1 := &projP1xP1{}
	tmp2 := &projP2{}

	// Accumulate the odd components first
	v := NewIdentityPoint()
	for i := 1; i < 64; i += 2 {
		pointTable[i/2].SelectInto(multiple, digits[i])
		tmp1.AddAffine(v, multiple)
		v.fromP1xP1(tmp1)
	}

	// Multiply by 16
	tmp2.FromP3(v)       // tmp2 =    v in P2 coords
	tmp1.Double(tmp2)    // tmp1 =  2*v in P1xP1 coords
	tmp2.FromP1xP1(tmp1) // tmp2 =  2*v in P2 coords
	tmp1.Double(tmp2)    // tmp1 =  4*v in P1xP1 coords
	tmp2.FromP1xP1(tmp1) // tmp2 =  4*v in P2 coords
	tmp1.Double(tmp2)    // tmp1 =  8*v in P1xP1 coords
	tmp2.FromP1xP1(tmp1) // tmp2 =  8*v in P2 coords
	tmp1.Double(tmp2)    // tmp1 = 16*v in P1xP1 coords
	v.fromP1xP1(tmp1)    // now v = 16*(odd components)

	// Accumulate the even components
	for i := 0; i < 64; i += 2 {
		pointTable[i/2].SelectInto(multiple, digits[i])
		tmp1.AddAffine(v, multiple)
		v.fromP1xP1(tmp1)
	}

	return v
}

// ScalarMultPrecomp sets v = x * P and returns v.
//
// The scalar multiplication is done in constant time.
// This is slightly slower than ScalarBaseMult, but much faster than ScalarMult
// for points which are used many times
func (v *Point) ScalarMultPrecomp(x *Scalar, p *Point) *Point {
	pointTable := PrecompPoint(p)

	v.Set(pointTable.ScalarMult(x))

	return v
}
