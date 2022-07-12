// Copyright (c) 2019 The Go Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package edwards25519

import (
	"sync"
)

// basepointTable is a set of 32 affineLookupTables, where table i is generated
// from 256i * basepoint. It is precomputed the first time it's used.
func basepointTable() *PrecomputedPoint {
	basepointTablePrecomp.initOnce.Do(func() {
		basepointTablePrecomp.table = PrecompPoint(NewGeneratorPoint())
	})
	return basepointTablePrecomp.table
}

var basepointTablePrecomp struct {
	table    *PrecomputedPoint
	initOnce sync.Once
}

// ScalarBaseMult sets v = x * B, where B is the canonical generator, and
// returns v.
//
// The scalar multiplication is done in constant time.
func (v *Point) ScalarBaseMult(x *Scalar) *Point {
	basepointTable := basepointTable()
	v.Set(basepointTable.ScalarMult(x))
	return v
}

// ScalarMult sets v = x * q, and returns v.
//
// The scalar multiplication is done in constant time.
func (v *Point) ScalarMult(x *Scalar, q *Point) *Point {
	checkInitialized(q)

	var table projLookupTable
	table.FromP3(q)

	// Write x = sum(x_i * 16^i)
	// so  x*Q = sum( Q*x_i*16^i )
	//         = Q*x_0 + 16*(Q*x_1 + 16*( ... + Q*x_63) ... )
	//           <------compute inside out---------
	//
	// We use the lookup table to get the x_i*Q values
	// and do four doublings to compute 16*Q
	digits := x.signedRadix16()

	// Unwrap first loop iteration to save computing 16*identity
	multiple := &projCached{}
	tmp1 := &projP1xP1{}
	tmp2 := &projP2{}
	table.SelectInto(multiple, digits[63])

	v.Set(NewIdentityPoint())
	tmp1.Add(v, multiple) // tmp1 = x_63*Q in P1xP1 coords
	for i := 62; i >= 0; i-- {
		tmp2.FromP1xP1(tmp1) // tmp2 =    (prev) in P2 coords
		tmp1.Double(tmp2)    // tmp1 =  2*(prev) in P1xP1 coords
		tmp2.FromP1xP1(tmp1) // tmp2 =  2*(prev) in P2 coords
		tmp1.Double(tmp2)    // tmp1 =  4*(prev) in P1xP1 coords
		tmp2.FromP1xP1(tmp1) // tmp2 =  4*(prev) in P2 coords
		tmp1.Double(tmp2)    // tmp1 =  8*(prev) in P1xP1 coords
		tmp2.FromP1xP1(tmp1) // tmp2 =  8*(prev) in P2 coords
		tmp1.Double(tmp2)    // tmp1 = 16*(prev) in P1xP1 coords
		v.fromP1xP1(tmp1)    //    v = 16*(prev) in P3 coords
		table.SelectInto(multiple, digits[i])
		tmp1.Add(v, multiple) // tmp1 = x_i*Q + 16*(prev) in P1xP1 coords
	}
	v.fromP1xP1(tmp1)
	return v
}

// basepointNafTable is the nafLookupTable8 for the basepoint.
// It is precomputed the first time it's used.
func basepointNafTable() *nafLookupTable8 {
	basepointNafTablePrecomp.initOnce.Do(func() {
		basepointNafTablePrecomp.table.FromP3(NewGeneratorPoint())
	})
	return &basepointNafTablePrecomp.table
}

var basepointNafTablePrecomp struct {
	table    nafLookupTable8
	initOnce sync.Once
}

// VarTimeDoubleScalarBaseMult sets v = a * A + b * B, where B is the canonical
// generator, and returns v.
//
// Execution time depends on the inputs.
func (v *Point) VarTimeDoubleScalarBaseMult(a *Scalar, A *Point, b *Scalar) *Point {
	checkInitialized(A)

	// Similarly to the single variable-base approach, we compute
	// digits and use them with a lookup table.  However, because
	// we are allowed to do variable-time operations, we don't
	// need constant-time lookups or constant-time digit
	// computations.
	//
	// So we use a non-adjacent form of some width w instead of
	// radix 16.  This is like a binary representation (one digit
	// for each binary place) but we allow the digits to grow in
	// magnitude up to 2^{w-1} so that the nonzero digits are as
	// sparse as possible.  Intuitively, this "condenses" the
	// "mass" of the scalar onto sparse coefficients (meaning
	// fewer additions).

	basepointNafTable := basepointNafTable()
	var aTable nafLookupTable5
	aTable.FromP3(A)
	// Because the basepoint is fixed, we can use a wider NAF
	// corresponding to a bigger table.
	aNaf := a.nonAdjacentForm(5)
	bNaf := b.nonAdjacentForm(8)

	// Find the first nonzero coefficient.
	i := 255
	for j := i; j >= 0; j-- {
		if aNaf[j] != 0 || bNaf[j] != 0 {
			break
		}
	}

	multA := &projCached{}
	multB := &affineCached{}
	tmp1 := &projP1xP1{}
	tmp2 := &projP2{}
	tmp2.Zero()

	// Move from high to low bits, doubling the accumulator
	// at each iteration and checking whether there is a nonzero
	// coefficient to look up a multiple of.
	for ; i >= 0; i-- {
		tmp1.Double(tmp2)

		// Only update v if we have a nonzero coeff to add in.
		if aNaf[i] > 0 {
			v.fromP1xP1(tmp1)
			aTable.SelectInto(multA, aNaf[i])
			tmp1.Add(v, multA)
		} else if aNaf[i] < 0 {
			v.fromP1xP1(tmp1)
			aTable.SelectInto(multA, -aNaf[i])
			tmp1.Sub(v, multA)
		}

		if bNaf[i] > 0 {
			v.fromP1xP1(tmp1)
			basepointNafTable.SelectInto(multB, bNaf[i])
			tmp1.AddAffine(v, multB)
		} else if bNaf[i] < 0 {
			v.fromP1xP1(tmp1)
			basepointNafTable.SelectInto(multB, -bNaf[i])
			tmp1.SubAffine(v, multB)
		}

		tmp2.FromP1xP1(tmp1)
	}

	v.fromP2(tmp2)
	return v
}
