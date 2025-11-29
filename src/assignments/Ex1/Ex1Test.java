package assignments.Ex1;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}
	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
        double[] po_a = {2, 1, -0.7, -0.02, 0.02};
        double[] po_b = {6, 0.1, -0.2};
        double x1 = Ex1.sameValue(po_a, po_b, -10, -5, Ex1.EPS);
        double a1 = Ex1.area(po_a, po_b, x1, 6, 8);
        double area = 58.5658;
        assertEquals(a1, area, Ex1.EPS);
    }
        @Test
        void testF_SimplePolynomial() {
            //cheak a simple polynom 2x+2
            //  x=0: f(0) =2+2*0=2
            //  x=1: f(1) = 2 + 2*1 = 4
            //  x=2: f(2) = 2 + 2*2 = 6
            double[] poly = {2, 2};
            assertEquals(2, Ex1.f(poly, 0), Ex1.EPS);
            assertEquals(4, Ex1.f(poly, 1), Ex1.EPS);
            assertEquals(6, Ex1.f(poly, 2), Ex1.EPS);
        }

        @Test
        void testF_QuadraticPolynomial() {
            // chek quarter poly י: 1 + 2x + 3x^2
            //  x=2: f(2) = 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17
            double[] poly = {1, 2, 3};
            assertEquals(17, Ex1.f(poly, 2), Ex1.EPS);
        }

        @Test
        void testF_NegativeValues() {
            //check the polinom has nagtive value: -3 + 0.61x + 0.2x^2
            double[] poly = {-3, 0.61, 0.2};
            double result = Ex1.f(poly, 1);
            // f(1) = -3 + 0.61 + 0.2 = -2.19
            assertEquals(-2.19, result, Ex1.EPS);
        }

        @Test
        void testF_ZeroPolynomial() {
        // cheak the output is 0
            assertEquals(0, Ex1.f(Ex1.ZERO, 5), Ex1.EPS);
            assertEquals(0, Ex1.f(Ex1.ZERO, -10), Ex1.EPS);
        }

        // ========== Test add ==========

        @Test
        void testAdd_BasicAddition() {
            // check simple addition: (2 + 2x) + (-3 + 0.61x + 0.2x^2)
            //  -1 + 2.61x + 0.2x^2
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double[] result = Ex1.add(p1, p2);

            assertEquals(-1, result[0], Ex1.EPS);
            assertEquals(2.61, result[1], Ex1.EPS);
            assertEquals(0.2, result[2], Ex1.EPS);
        }

        @Test
        void testAdd_Commutative() {
            // check commutative: p1 + p2 = p2 + p1
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double[] sum1 = Ex1.add(p1, p2);
            double[] sum2 = Ex1.add(p2, p1);

            assertTrue(Ex1.equals(sum1, sum2));
        }

        @Test
        void testAdd_WithZero() {
            //  p + 0 = p
            double[] p = {2, 2};
            double[] result = Ex1.add(p, Ex1.ZERO);
            assertTrue(Ex1.equals(p, result));
        }

        @Test
        void testAdd_InverseProperty() {
            // cheak p + (-p) = 0
            double[] p = {2, 2};
            double[] minusOne = {-1};
            double[] negativeP = Ex1.mul(p, minusOne);
            double[] sum = Ex1.add(p, negativeP);

            assertTrue(Ex1.equals(sum, Ex1.ZERO));
        }

        @Test
        void testAdd_DifferentLength(){
        // cheak the sum of polinom who has different length
            double[] p1 = {1}; // 1
            double[] p2 = {2, 3, 4, 5}; // 2 + 3x + 4x^2 + 5x^3
            double[] result = Ex1.add(p1, p2);

            assertEquals(3, result[0], Ex1.EPS); // 1+2
            assertEquals(3, result[1], Ex1.EPS);
            assertEquals(4, result[2], Ex1.EPS);
            assertEquals(5, result[3], Ex1.EPS);
        }

        @Test
        void testAdd_BothZero() {
            // cheak if add to two polynom =0
            double[] result = Ex1.add(Ex1.ZERO, Ex1.ZERO);
            assertTrue(Ex1.equals(result, Ex1.ZERO));
        }

        @Test
        void testAdd_Associative() {
            //  (p1 + p2) + p3 = p1 + (p2 + p3)
            double[] p1 = {1, 2};
            double[] p2 = {3, 4};
            double[] p3 = {5, 6};

            double[] left = Ex1.add(Ex1.add(p1, p2), p3);
            double[] right = Ex1.add(p1, Ex1.add(p2, p3));

            assertTrue(Ex1.equals(left, right));
        }

        @Test
        void testAdd_NegativeCoefficients() {
            //cheak to add the polynom who has nagtive value
            double[] p1 = {-5, -3, -1};
            double[] p2 = {5, 3, 1};
            double[] result = Ex1.add(p1, p2);

            assertTrue(Ex1.equals(result, Ex1.ZERO));
        }

        @Test
        void testAdd_VerifyWithF() {
            //  x: f(p1+p2, x) = f(p1, x) + f(p2, x)
            double[] p1 = {1, 2, 3};
            double[] p2 = {4, 5};
            double[] sum = Ex1.add(p1, p2);
            double x = Math.PI;
            double expected = Ex1.f(p1, x) + Ex1.f(p2, x);
            double actual = Ex1.f(sum, x);

            assertEquals(expected, actual ,Ex1.EPS);
        }

// test for mul//
        @Test
        void testMul_BasicMultiplication() {
            //  (2 + 2x) * (-3 + 0.61x + 0.2x^2)
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double[] result = Ex1.mul(p1, p2);


            assertEquals(-6, result[0], Ex1.EPS);  // 2*(-3)
            assertEquals(-4.78, result[1], Ex1.EPS); // 2*0.61 + 2*(-3)
        }

        @Test
        void testMul_Commutative() {
            //: p1 * p2 = p2 * p1
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double[] prod1 = Ex1.mul(p1, p2);
            double[] prod2 = Ex1.mul(p2, p1);

            assertTrue(Ex1.equals(prod1, prod2));
        }

        @Test
        void testMul_WithZero() {
            //  0: p * 0 = 0
            double[] p = {2, 2};
            double[] result = Ex1.mul(p, Ex1.ZERO);
            assertTrue(Ex1.equals(result, Ex1.ZERO));
        }

        @Test
        void testMul_VerifyWithF() {
            //  x: f(p1*p2, x) = f(p1, x) * f(p2, x)
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double[] product = Ex1.mul(p1, p2);

            double[] testPoints = {0, 1, 2, 3, -1, -5};
            for (double x : testPoints) {
                double expected = Ex1.f(p1, x) * Ex1.f(p2, x);
                double actual = Ex1.f(product, x);
                assertEquals(expected, actual, Ex1.EPS);
            }
        }

        @Test
        void testMul_ByConstant() {
            // : (1 + 2x + 3x^2) * 5
            double[] poly = {1, 2, 3};
            double[] constant = {5};
            double[] result = Ex1.mul(poly, constant);

            assertEquals(5, result[0], Ex1.EPS);
            assertEquals(10, result[1], Ex1.EPS);
            assertEquals(15, result[2], Ex1.EPS);
        }

        @Test
        void testMul_TwoLinear() {
            // multiple two polynom (1 + x) * (2 + x) = 2 + 3x + x^2
            double[] p1 = {1, 1};
            double[] p2 = {2, 1};
            double[] result = Ex1.mul(p1, p2);

            assertEquals(2, result[0], Ex1.EPS);
            assertEquals(3, result[1], Ex1.EPS);
            assertEquals(1, result[2], Ex1.EPS);
        }

        @Test
        void testMul_WithOne() {
            // cheak mul by 1: p * 1 = p
            double[] poly = {1, 2, 3};
            double[] one = {1};
            double[] result = Ex1.mul(poly, one);

            assertTrue(Ex1.equals(poly, result));
        }

        @Test
        void testMul_ResultLength() {
            // length(p1*p2) = length(p1) + length(p2) - 1
            double[] p1 = {1, 2, 3}; // degree 2
            double[] p2 = {4, 5};    // degree 1
            double[] result = Ex1.mul(p1, p2);

            assertEquals(4, result.length); // 3 + 2 - 1 = 4
        }

        @Test
        void testMul_NegativeCoefficients() {
        //multiple twopolynom who has negitve value
            double[] p1 = {-1, 2};  // -1 + 2x
            double[] p2 = {3, -1};  // 3 - x
            double[] result = Ex1.mul(p1, p2);

            assertEquals(-3, result[0], Ex1.EPS);  // -1*3
            assertEquals(7, result[1], Ex1.EPS);   // -1*(-1) + 2*3
            assertEquals(-2, result[2], Ex1.EPS);  // 2*(-1)
        }

        @Test
        void testMul_Associative() {
            // : (p1 * p2) * p3 = p1 * (p2 * p3)
            double[] p1 = {1, 1};
            double[] p2 = {1, 1};
            double[] p3 = {1, 1};

            double[] left = Ex1.mul(Ex1.mul(p1, p2), p3);
            double[] right = Ex1.mul(p1, Ex1.mul(p2, p3));

            assertTrue(Ex1.equals(left, right));
        }

        // ========== test derivative ==========

        @Test
        void testDerivative_Quadratic() {
            //  1 + 2x + 3x^2 היא 2 + 6x
            double[] poly = {1, 2, 3};
            double[] derivative = Ex1.derivative(poly);

            assertEquals(2, derivative.length);
            assertEquals(2, derivative[0], Ex1.EPS);
            assertEquals(6, derivative[1], Ex1.EPS);
        }

        @Test
        void testDerivative_Linear() {
            //  5 + 3x +3
            double[] poly = {5, 3};
            double[] derivative = Ex1.derivative(poly);

            assertEquals(1, derivative.length);
            assertEquals(3, derivative[0], Ex1.EPS);
        }

        @Test
        void testDerivative_Constant() {
            //0 derivative
            double[] poly = {7};
            double[] derivative = Ex1.derivative(poly);
            assertTrue(Ex1.equals(derivative, Ex1.ZERO));
        }

        @Test
        void testDerivative_RepeatedDerivatives() {
            // 0
            double[] poly = {1, 2, 3}; // 1 + 2x + 3x^2
            double[] d1 = Ex1.derivative(poly);  // 2 + 6x
            double[] d2 = Ex1.derivative(d1);    // 6
            double[] d3 = Ex1.derivative(d2);    // 0
            double[] d4 = Ex1.derivative(d3);    // 0

            assertTrue(Ex1.equals(d3, Ex1.ZERO));
            assertTrue(Ex1.equals(d4, Ex1.ZERO));
        }

        @Test
        void testDerivative_HighDegree() {
            // : x^5
            double[] poly = {0, 0, 0, 0, 0, 1}; // x^5
            double[] derivative = Ex1.derivative(poly);

            assertEquals(5, derivative.length);
            assertEquals(0, derivative[0], Ex1.EPS);
            assertEquals(0, derivative[1], Ex1.EPS);
            assertEquals(0, derivative[2], Ex1.EPS);
            assertEquals(0, derivative[3], Ex1.EPS);
            assertEquals(5, derivative[4], Ex1.EPS); // 5x^4
        }

        @Test
        void testDerivative_WithNegativeCoefficients() {
            // : -3 + 2x - x^2
            double[] poly = {-3, 2, -1};
            double[] derivative = Ex1.derivative(poly);

            assertEquals(2, derivative[0], Ex1.EPS);
            assertEquals(-2, derivative[1], Ex1.EPS);
        }

        @Test
        void testDerivative_ZeroPolynomial() {
            // derivative 0 is 0
            double[] derivative = Ex1.derivative(Ex1.ZERO);
            assertTrue(Ex1.equals(derivative, Ex1.ZERO));
        }

        @Test
        void testDerivative_SumRule() {
            //: (f+g)' = f' + g'
            double[] p1 = {1, 2, 3};
            double[] p2 = {4, 5, 6};

            double[] sumDeriv = Ex1.derivative(Ex1.add(p1, p2));
            double[] derivSum = Ex1.add(Ex1.derivative(p1), Ex1.derivative(p2));

            assertTrue(Ex1.equals(sumDeriv, derivSum));
        }

        @Test
        void testDerivative_ProductRule_Simple() {
            //  x * x = x^2
            // (x*x)' = 2x,  : x'*x + x*x' = 1*x + x*1 = 2x
            double[] x = {0, 1}; // x
            double[] x_squared = Ex1.mul(x, x); // x^2
            double[] deriv = Ex1.derivative(x_squared); // 2x

            assertEquals(2, deriv[0], Ex1.EPS);
            assertEquals(2, deriv[1], Ex1.EPS);
        }

        // ==========  Test equals ==========

        @Test
        void testEquals_SamePolynomial() {
            // polynom equal to itself
            double[] p = {1, 2, 3};
            assertTrue(Ex1.equals(p, p));
        }

        @Test
        void testEquals_WithTrailingZeros() {
            // {1, 2, 0, 0} שווה ל-{1, 2}
            double[] p1 = {1, 2, 0, 0};
            double[] p2 = {1, 2};
            assertTrue(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_WithinEpsilon() {
            // epsilon
            double[] p1 = {1.0};
            double[] p2 = {1.0 + Ex1.EPS / 2};
            assertTrue(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_OutsideEpsilon(){
            double[] p1 = {1.0};
            double[] p2 = {1.0 + Ex1.EPS * 2};
            assertFalse(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_BothZero() {
            // two polinom =0
            assertTrue(Ex1.equals(Ex1.ZERO, Ex1.ZERO));
            assertTrue(Ex1.equals(Ex1.ZERO, new double[]{0}));
            assertTrue(Ex1.equals(new double[]{0, 0, 0}, Ex1.ZERO));
        }

        @Test
        void testEquals_DifferentDegrees() {
            //
            double[] p1 = {1, 2};     // 1 + 2x
            double[] p2 = {1, 2, 1};  // 1 + 2x + x^2
            assertFalse(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_Symmetric() {
            // : equals(p1,p2) = equals(p2,p1)
            double[] p1 = {1, 2, 3};
            double[] p2 = {1, 2, 3};
            assertTrue(Ex1.equals(p1, p2));
            assertTrue(Ex1.equals(p2, p1));
        }

        @Test
        void testEquals_NegativeCoefficients() {
            // two polynom  has negetive value
            double[] p1 = {-1, -2, -3};
            double[] p2 = {-1, -2, -3};
            assertTrue(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_SmallDifferences() {
            // cheak  the different between
            double[] p1 = {1.0000001, 2.0000001};
            double[] p2 = {1.0, 2.0};
            assertTrue(Ex1.equals(p1, p2));
        }

        @Test
        void testEquals_OneCoefficientDifferent() {
            //  only one value isnt same
            double[] p1 = {1, 2, 3};
            double[] p2 = {1, 2, 4};
            assertFalse(Ex1.equals(p1, p2));
        }

        // =========  test poly ==========

        @Test
        void testPoly_SimplePolynomial() {
            //  2 + 2x
            double[] poly = {2, 2};
            String result = Ex1.poly(poly);
            assertTrue(result.contains("2.0x") || result.contains("2x"));
            assertTrue(result.contains("2.0") || result.contains("+2"));
        }

        @Test
        void testPoly_WithNegativeCoefficients() {
            //cheak negtive
            double[] poly = {-1.1, 2.3, 3.1};
            String result = Ex1.poly(poly);
            assertTrue(result.contains("3.1"));
            assertTrue(result.contains("-1.1"));
        }

        @Test
        void testPoly_SkipsZeroCoefficients() {
            // if 0 is value
            double[] poly = {1, 0, 3};
            String result = Ex1.poly(poly);
            // we dont need write 0x
            assertFalse(result.contains("0x") || result.contains("0.0x"));
        }

        @Test
        void testPoly_ConstantOnly() {
            // בודק פולינום קבוע
            double[] poly = {5.5};
            String result = Ex1.poly(poly);
            assertTrue(result.contains("5.5"));
            assertFalse(result.contains("x"));
        }

        @Test
        void testPoly_LinearOnly() {
            // בודק פולינום לינארי: 0 + 3x
            double[] poly = {0, 3};
            String result = Ex1.poly(poly);
            assertTrue(result.contains("3") && result.contains("x"));
        }

        @Test
        void testPoly_HighDegree() {
            // בודק פולינום ממעלה גבוהה
            double[] poly = {1, 0, 0, 0, 2}; // 1 + 2x^4
            String result = Ex1.poly(poly);
            assertTrue(result.contains("x^4"));
            assertTrue(result.contains("1"));
        }

        @Test
        void testPoly_AllNegative() {
            // בודק פולינום עם כל המקדמים שליליים
            double[] poly = {-1, -2, -3};
            String result = Ex1.poly(poly);
            assertTrue(result.contains("-1"));
            assertTrue(result.contains("-2"));
            assertTrue(result.contains("-3"));
        }

        @Test
        void testPoly_ZeroPolynomial() {
            // בודק פולינום אפס
            String result = Ex1.poly(Ex1.ZERO);
            assertTrue(result.equals("0") || result.isEmpty());
        }

        @Test
        void testPoly_EmptyArray() {
            // בודק מערך ריק
            double[] empty = {};
            String result = Ex1.poly(empty);
            assertEquals("0", result);
        }

        // ========== טסטים לפונקציה getPolynomFromString ==========

        @Test
        void testGetPolynomFromString_Simple() {
            // בודק פרסור של מחרוזת פשוטה
            String input = "3.1x^2 +2.3x -1.1";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(-1.1, result[0], Ex1.EPS);
            assertEquals(2.3, result[1], Ex1.EPS);
            assertEquals(3.1, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_RoundTrip() {
            // בודק שהמרה למחרוזת וחזרה נותנת את אותו פולינום
            double[] original = {-1.1, 2.3, 3.1};
            String str = Ex1.poly(original);
            double[] parsed = Ex1.getPolynomFromString(str);

            assertTrue(Ex1.equals(original, parsed));
        }

        @Test
        void testGetPolynomFromString_WithSpaces() {
            // בודק פרסור עם רווחים
            String input = "2x^2 + 3x - 1";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(-1, result[0], Ex1.EPS);
            assertEquals(3, result[1], Ex1.EPS);
            assertEquals(2, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_ConstantOnly() {
            // בודק פרסור של קבוע בלבד
            String input = "5";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(1, result.length);
            assertEquals(5, result[0], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_NoSpaces() {
            // בודק פרסור ללא רווחים
            String input = "3x^2+2x-1";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(-1, result[0], Ex1.EPS);
            assertEquals(2, result[1], Ex1.EPS);
            assertEquals(3, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_NegativeFirst() {
            // בודק פרסור כשהמקדם הראשון שלילי
            String input = "-x^2 +x +1";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(1, result[0], Ex1.EPS);
            assertEquals(1, result[1], Ex1.EPS);
            assertEquals(-1, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_ImplicitOne() {
            // בודק פרסור כשהמקדם הוא 1 (לא רשום)
            String input = "x^2 +x";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(0, result[0], Ex1.EPS);
            assertEquals(1, result[1], Ex1.EPS);
            assertEquals(1, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_ImplicitMinusOne() {
            // בודק פרסור כשהמקדם הוא -1 (רק מינוס)
            String input = "-x^2 -x";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(0, result[0], Ex1.EPS);
            assertEquals(-1, result[1], Ex1.EPS);
            assertEquals(-1, result[2], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_Zero() {
            // בודק פרסור של אפס
            String input = "0";
            double[] result = Ex1.getPolynomFromString(input);
            assertTrue(Ex1.equals(result, Ex1.ZERO));
        }

        @Test
        void testGetPolynomFromString_EmptyString() {
            // בודק פרסור של מחרוזת ריקה
            String input = "";
            double[] result = Ex1.getPolynomFromString(input);
            assertTrue(Ex1.equals(result, Ex1.ZERO));
        }

        @Test
        void testGetPolynomFromString_HighDegree() {
            // בודק פרסור של מעלה גבוהה
            String input = "2x^5 -3x^3 +1";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(6, result.length);
            assertEquals(1, result[0], Ex1.EPS);
            assertEquals(0, result[1], Ex1.EPS);
            assertEquals(0, result[2], Ex1.EPS);
            assertEquals(-3, result[3], Ex1.EPS);
            assertEquals(0, result[4], Ex1.EPS);
            assertEquals(2, result[5], Ex1.EPS);
        }

        @Test
        void testGetPolynomFromString_DecimalCoefficients() {
            // בודק פרסור עם מקדמים עשרוניים
            String input = "1.5x^2 +0.5x -2.3";
            double[] result = Ex1.getPolynomFromString(input);

            assertEquals(-2.3, result[0], Ex1.EPS);
            assertEquals(0.5, result[1], Ex1.EPS);
            assertEquals(1.5, result[2], Ex1.EPS);
        }

        // ========== טסטים לפונקציה root_rec ==========

        @Test
        void testRootRec_SimpleRoot() {
            // בודק מציאת שורש של x - 2 = 0 (שורש ב-x=2)
            double[] poly = {-2, 1}; // x - 2
            double root = Ex1.root_rec(poly, 0, 5, Ex1.EPS);
            assertEquals(2, root, Ex1.EPS);
        }

        @Test
        void testRootRec_QuadraticRoot() {
            // בודק מציאת שורש של פולינום ריבועי
            double[] poly = {-2, 0, 1}; // x^2 - 2 (שורש ב-√2)
            double root = Ex1.root_rec(poly, 0, 2, Ex1.EPS);
            assertEquals(Math.sqrt(2), root, Ex1.EPS);
        }

        @Test
        void testRootRec_NegativeRoot() {
            // בודק מציאת שורש שלילי: x + 3 = 0 (שורש ב-x=-3)
            double[] poly = {3, 1}; // x + 3
            double root = Ex1.root_rec(poly, -5, 0, Ex1.EPS);
            assertEquals(-3, root, Ex1.EPS);
        }

        @Test
        void testRootRec_AtZero() {
            // בודק מציאת שורש ב-0: x = 0
            double[] poly = {0, 1}; // x
            double root = Ex1.root_rec(poly, -1, 1, Ex1.EPS);
            assertEquals(0, root, Ex1.EPS);
        }

        @Test
        void testRootRec_CubicRoot() {
            // בודק מציאת שורש של פולינום מעוקב: x^3 - 8 = 0 (שורש ב-2)
            double[] poly = {-8, 0, 0, 1}; // x^3 - 8
            double root = Ex1.root_rec(poly, 0, 5, Ex1.EPS);
            assertEquals(2, root, 0.01);
        }

        @Test
        void testRootRec_FractionRoot() {
            // בודק מציאת שורש שבר: 2x - 1 = 0 (שורש ב-0.5)
            double[] poly = {-1, 2}; // 2x - 1
            double root = Ex1.root_rec(poly, 0, 1, Ex1.EPS);
            assertEquals(0.5, root, Ex1.EPS);
        }

        // ========== טסטים לפונקציה sameValue ==========

        @Test
        void testSameValue_Symmetric() {
            // בודק שהפונקציה סימטרית: sameValue(p1,p2) = sameValue(p2,p1)
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double x1 = -4, x2 = 0;

            double result1 = Ex1.sameValue(p1, p2, x1, x2, Ex1.EPS);
            double result2 = Ex1.sameValue(p2, p1, x1, x2, Ex1.EPS);

            assertEquals(result1, result2, Ex1.EPS);
        }

        @Test
        void testSameValue_TwoLines() {
            // בודק נקודת חיתוך של שני קווים: y=x ו-y=2-x
            // נפגשים ב-x=1
            double[] p1 = {0, 1};   // y = x
            double[] p2 = {2, -1};  // y = 2 - x
            double intersection = Ex1.sameValue(p1, p2, 0, 3, Ex1.EPS);

            assertEquals(1, intersection, Ex1.EPS);
        }

        @Test
        void testSameValue_ParabolaAndLine() {
            // בודק חיתוך של פרבולה וקו: y=x^2 ו-y=x
            // נפגשים ב-x=0 או x=1
            double[] p1 = {0, 0, 1};  // y = x^2
            double[] p2 = {0, 1};     // y = x
            double intersection = Ex1.sameValue(p1, p2, 0.5, 2, Ex1.EPS);

            assertEquals(1, intersection, Ex1.EPS);
        }

        @Test
        void testSameValue_SamePolynomial() {
            // בודק שני פולינומים זהים - כל נקודה היא נקודת "חיתוך"
            double[] p = {1, 2, 3};
            double result = Ex1.sameValue(p, p, 0, 5, Ex1.EPS);

            // כל נקודה בטווח תהיה תקינה
            assertTrue(result >= 0 && result <= 5);
        }

        @Test
        void testSameValue_NegativeRange() {
            // בודק חיתוך בטווח שלילי
            double[] p1 = {0, 1};   // y = x
            double[] p2 = {0, -1};  // y = -x
            double intersection = Ex1.sameValue(p1, p2, -2, 2, Ex1.EPS);

            assertEquals(0, intersection, Ex1.EPS);
        }

        // ========== טסטים לפונקציה length ==========

        @Test
        void testLength_StraightLine() {
            // בודק אורך של קו ישר y = x מ-0 ל-1
            // אורך צריך להיות √2 ≈ 1.414
            double[] poly = {0, 1}; // y = x
            double length = Ex1.length(poly, 0, 1, 100);
            assertEquals(Math.sqrt(2), length, 0.01);
        }

        @Test
        void testLength_Constant() {
            // בודק אורך של פונקציה קבועה (קו אופקי)
            double[] poly = {5}; // y = 5
            double length = Ex1.length(poly, 0, 10, 100);
            // אורך צריך להיות 10 (רק בכיוון x)
            assertEquals(10, length, Ex1.EPS);
        }

        @Test
        void testLength_VerticalChange() {
            // בודק אורך של פונקציה עם שינוי אנכי בלבד נטו
            // למשל 2x מ-0 ל-5: אורך צריך להיות גדול מ-10
            double[] poly = {0, 2}; // y = 2x
            double length = Ex1.length(poly, 0, 5, 100);

            // אורך הפונקציה: √(5² + 10²) = √125 ≈ 11.18
            assertEquals(Math.sqrt(125), length, 0.1);
        }

        @Test
        void testLength_Parabola() {
            // בודק אורך של פרבולה y = x^2 מ-0 ל-1
            double[] poly = {0, 0, 1}; // y = x^2
            double length = Ex1.length(poly, 0, 1, 1000);

            // האורך המדויק הוא בערך 1.478
            assertEquals(1.478, length, 0.01);
        }

        @Test
        void testLength_NegativeRange() {
            // בודק אורך בטווח שלילי
            double[] poly = {0, 1}; // y = x
            double length = Ex1.length(poly, -2, -1, 100);

            assertEquals(Math.sqrt(2), length, 0.01);
        }

        @Test
        void testLength_MoreSegmentsBetterAccuracy() {
            // בודק שיותר קטעים נותן דיוק טוב יותר
            double[] poly = {0, 0, 1}; // y = x^2

            double length10 = Ex1.length(poly, 0, 1, 10);
            double length100 = Ex1.length(poly, 0, 1, 100);
            double length1000 = Ex1.length(poly, 0, 1, 1000);

            // כל חישוב עם יותר קטעים צריך להיות יותר מדויק
            assertTrue(Math.abs(length1000 - 1.478) < Math.abs(length100 - 1.478));
            assertTrue(Math.abs(length100 - 1.478) < Math.abs(length10 - 1.478));
        }

        @Test
        void testLength_ZeroLength() {
            // בודק אורך כאשר x1 = x2
            double[] poly = {1, 2, 3};
            double length = Ex1.length(poly, 5, 5, 100);

            assertEquals(0, length, Ex1.EPS);
        }

        // ========== טסטים לפונקציה area ==========

        @Test
        void testArea_Symmetric() {
            // בודק שחישוב שטח סימטרי: area(p1,p2) = area(p2,p1)
            double[] p1 = {2, 2};
            double[] p2 = {-3, 0.61, 0.2};
            double x1 = -4, x2 = 0;

            double area1 = Ex1.area(p1, p2, x1, x2, 100);
            double area2 = Ex1.area(p2, p1, x1, x2, 100);

            assertEquals(area1, area2, Ex1.EPS);
        }

        @Test
        void testArea_SimpleTriangle() {
            // בודק שטח בין y=0 ל-y=x בטווח [-1,2]
            // זה משולש עם בסיס 3 וגובה 2, שטח = 2.5
            double[] p1 = Ex1.ZERO;
            double[] p2 = {0, 1};

            double area = Ex1.area(p1, p2, -1, 2, 100);
            assertEquals(2.5, area, Ex1.EPS);
        }

        @Test
        void testArea_WithIntersections() {
            // בודק שטח כאשר הפונקציות נחתכות בתוך הטווח
            double[] p1 = {2, 1, -0.7, -0.02, 0.02};
            double[] p2 = {6, 0.1, -0.2};

            double x1 = Ex1.sameValue(p1, p2, -10, -5, Ex1.EPS);
            double area = Ex1.area(p1, p2, x1, 6, 8);

            assertEquals(58.5658, area, Ex1.EPS);
        }

        @Test
        void testArea_Rectangle() {
            // בודק שטח מלבן: בין y=0 ל-y=2 מ-x=0 ל-x=5
            // שטח = 2 * 5 = 10
            double[] p1 = Ex1.ZERO;
            double[] p2 = {2};  // y = 2

            double area = Ex1.area(p1, p2, 0, 5, 100);
            assertEquals(10, area, Ex1.EPS);
        }

        @Test
        void testArea_TwoLines() {
            // בודק שטח בין שני קווים: y=x ו-y=2x מ-0 ל-2
            // השטח בניהם: ∫(2x-x)dx = ∫x dx = x²/2 מ-0 ל-2 = 2
            double[] p1 = {0, 1};   // y = x
            double[] p2 = {0, 2};   // y = 2x

            double area = Ex1.area(p1, p2, 0, 2, 100);
            assertEquals(2, area, 0.01);
        }

        @Test
        void testArea_ParabolaAndLine() {
            // בודק שטח בין y=x^2 ו-y=x מ-0 ל-1
            // ∫(x-x²)dx = x²/2 - x³/3 מ-0 ל-1 = 1/2 - 1/3 = 1/6 ≈ 0.167
            double[] p1 = {0, 0, 1};  // y = x^2
            double[] p2 = {0, 1};     // y = x

            double area = Ex1.area(p1, p2, 0, 1, 100);
            assertEquals(1.0/6.0, area, 0.01);
        }

        @Test
        void testArea_NegativeRange() {
            // בודק שטח בטווח שלילי
            double[] p1 = Ex1.ZERO;
            double[] p2 = {0, 1};  // y = x

            double area = Ex1.area(p1, p2, -2, -1, 100);
            // |∫x dx| מ--2 ל--1 = |x²/2| = |1/2 - 2| = 1.5
            assertEquals(1.5, area, 0.01);
        }

        @Test
        void testArea_SamePolynomial() {
            // בודק שטח בין פולינום לעצמו - צריך להיות 0
            double[] p = {1, 2, 3};
            double area = Ex1.area(p, p, 0, 5, 100);

            assertEquals(0, area, Ex1.EPS);
        }

        @Test
        void testArea_MoreTrapezoidsBetterAccuracy() {
            // בודק שיותר טרפזים נותן דיוק טוב יותר
            double[] p1 = {0, 0, 1};
            double[] p2 = {0, 1};

            double area10 = Ex1.area(p1, p2, 0, 1, 10);
            double area100 = Ex1.area(p1, p2, 0, 1, 100);
            double area1000 = Ex1.area(p1, p2, 0, 1, 1000);

            double expected = 1.0/6.0;
            assertTrue(Math.abs(area1000 - expected) < Math.abs(area100 - expected));
            assertTrue(Math.abs(area100 - expected) < Math.abs(area10 - expected));
        }

        @Test
        void testArea_CrossingFunctions() {
            // בודק שטח כאשר הפונקציות חוצות זו את זו
            // y=x ו-y=1-x נפגשות ב-x=0.5
            double[] p1 = {0, 1};   // y = x
            double[] p2 = {1, -1};  // y = 1-x

            double area = Ex1.area(p1, p2, 0, 1, 100);
            // שטח המשולש: 0.5 * 0.5 * 1 = 0.25
            assertEquals(0.25, area, 0.01);
        }

        // ========== טסטים לפונקציה PolynomFromPoints ==========

        @Test
        void testPolynomFromPoints_TwoPoints() {
            // בודק יצירת פולינום מ-2 נקודות (קו ישר)
            double[] xx = {0, 1};
            double[] yy = {2, 4};  // y = 2x + 2

            double[] poly = Ex1.PolynomFromPoints(xx, yy);

            assertNotNull(poly);
            assertEquals(2, poly.length);
            assertEquals(2, poly[0], Ex1.EPS);  // b = 2
            assertEquals(2, poly[1], Ex1.EPS);  // m = 2
        }

        @Test
        void testPolynomFromPoints_ThreePoints() {
            // בודק יצירת פולינום מ-3 נקודות (פרבולה)
            double[] xx = {0, 1, 2};
            double[] yy = {0, 1, 4};  // y = x^2

            double[] poly = Ex1.PolynomFromPoints(xx, yy);

            assertNotNull(poly);
            assertEquals(3, poly.length);
            // בודק שהפולינום עובר דרך הנקודות
            assertEquals(0, Ex1.f(poly, 0), 0.1);
            assertEquals(1, Ex1.f(poly, 1), 0.1);
            assertEquals(4, Ex1.f(poly, 2), 0.1);
        }

        @Test
        void testPolynomFromPoints_InvalidInput() {
            // בודק טיפול בקלט לא תקין
            assertNull(Ex1.PolynomFromPoints(null, null));
            assertNull(Ex1.PolynomFromPoints(new double[]{1}, new double[]{1, 2}));
            assertNull(Ex1.PolynomFromPoints(new double[]{1, 2, 3, 4}, new double[]{1, 2, 3, 4}));
        }
    }
