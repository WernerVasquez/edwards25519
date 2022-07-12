package edwards25519

import "testing"

func TestMultPrecompVsDalek(t *testing.T) {
	var p Point
	p.ScalarMultPrecomp(&dalekScalar, B)
	if dalekScalarBasepoint.Equal(&p) != 1 {
		t.Error("Scalar mul does not match dalek")
	}
	checkOnCurve(t, &p)
}

func TestMultPrecompVsDalekTwice(t *testing.T) {
	var p Point

	//fmt.Println(B.Bytes())

	p.ScalarMultPrecomp(&dalekScalar, B)
	if dalekScalarBasepoint.Equal(&p) != 1 {
		t.Error("Scalar mul does not match dalek")
	}

	//fmt.Println(B.Bytes())

	p.ScalarMultPrecomp(&dalekScalar, B)
	if dalekScalarBasepoint.Equal(&p) != 1 {
		t.Error("Scalar mul does not match dalek")
	}
}

func TestMultPrecompVsScalarMultTwoPoints(t *testing.T) {
	var p1, p2, p3, p4 Point
	p1.Set(dalekScalarBasepoint)
	p2.Set(dalekScalarBasepoint)

	p3.ScalarMultPrecomp(&dalekScalar, &p1)
	p4.ScalarMult(&dalekScalar, &p2)

	if p3.Equal(&p4) != 1 {
		t.Error("Scalar mul precomp does not match Scalar mul")
	}

	p1.Set(&p3)
	p2.Set(&p3)

	p3.ScalarMultPrecomp(&dalekScalar, &p1)
	p4.ScalarMult(&dalekScalar, &p2)

	if p3.Equal(&p4) != 1 {
		t.Error("Scalar mul precomp does not match Scalar mul")
	}
}

func TestPrecomputedPoint_ScalarMultVsDalek(t *testing.T) {

	p := PrecompPoint(B).ScalarMult(&dalekScalar)

	if dalekScalarBasepoint.Equal(p) != 1 {
		t.Error("Scalar mul does not match dalek")
	}
	checkOnCurve(t, p)
}

func TestScalarMultPrecompVsDalek(t *testing.T) {

	p1, precomputedPoint := new(Point).ScalarMultPrecomp(&dalekScalar, B)

	p2 := precomputedPoint.ScalarMult(&dalekScalar)

	if dalekScalarBasepoint.Equal(p1) != 1 {
		t.Error("ScalarMultPrecomp does not match dalek")
	}

	if dalekScalarBasepoint.Equal(p2) != 1 {
		t.Error("ScalarMult from returned precomputedPoint does not match dalek")
	}

	checkOnCurve(t, p1)
}

//benchmarks

func BenchmarkPrecompPoint(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PrecompPoint(NewGeneratorPoint())
	}
}

func BenchmarkPrecompPointInitialRunOnly(b *testing.B) {
	for i := 0; i < b.N; i++ {
		PrecompPoint(NewGeneratorPoint())

		//reset pointTablePrecomp so each run does full initialization
		pointTablePrecomp = make(map[[32]byte]PrecomputedPoint)
	}
}

func BenchmarkScalarMultPrecomp(b *testing.B) {
	var p Point

	for i := 0; i < b.N; i++ {
		p.ScalarMultPrecomp(&dalekScalar, B)
	}
}

func BenchmarkPrecomputedPoint_ScalarMult(b *testing.B) {

	t := PrecompPoint(B)

	for i := 0; i < b.N; i++ {
		t.ScalarMult(&dalekScalar)
	}
}
