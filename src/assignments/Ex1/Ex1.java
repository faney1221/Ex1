  package assignments.Ex1;

  import java.util.ArrayList;
  import java.util.Collections;

  /**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *a
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx
	 * @param yy
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans =null;
        if ( xx==null || yy==null ) {
            return null;
        }
        if (xx.length!=yy.length) { return null; }
        if (xx.length<2 || xx.length>3) { return null; }

		int lx = xx.length;
		int ly = yy.length;
        double x0=xx[0];
        double y0=yy[0];
        double x1=xx[1];
        double y1=yy[1];
		/** add you code below
		/////////////////// */
        if (xx.length==2) {
            double m=(y1-y0)/(x1-x0);
            double b = y1-m*x1;
            ans = new double[] {b,m};

        }
        if(xx.length==3) {
         double x2=xx[2];
         double y2=yy[2];
         double a=(y0*(x1-x2)+y1*(x2-x0)+y2*(x0-x1))/((x0-x1)*(x0-x2)*(x1-x2));
         double b =(y0-y1)/(x0-x1)-a*(x0-x1);
         double c=y0-a*x0*x0-b*y0;
         ans = new double[] {c,b,a};
        }

        return ans;


	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2.
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        /** add you code below
         */
        int maxd=Math.max(p1.length,p2.length);
        for (int i = 0; i <= maxd; i++) {
            double x =i;
            if(Math.abs( f(p1 ,x)-f(p2 ,x ))>EPS) {return false;}
        }
         /////////////////// */
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            /** add you code below
             */
for (int i=poly.length-1;i>=0;i--){
    if(poly[i]==0) continue;
    if (poly[i]>0 && !ans.isEmpty()) ans+= "+";

    if (poly[i]<0) ans+= "";

    if (i==0){
        ans += poly[i];}

     else if(i==1){ans += poly[i] +"x";
    }
     else {ans += poly[i]+"x^"+i;}
    }
}

        /////////////////// */
        return ans.trim();
		}

	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        /** add you code below
         */
        double[]diff=new double[Math.max(p1.length,p2.length)];
        for (int i=0;i<diff.length;i++) {
            double v1 = i < p1.length ? p1[i] : 0;
            double v2 = i < p2.length ? p2[i] : 0;
            diff[i] = v1 - v2;
        }
         /////////////////// */
		return root_rec(diff, x1, x2, eps);
	}
	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = 0;
        /** add you code below
         */
double step =(x2-x1)/numberOfSegments;
for (int i=0;i<numberOfSegments;i++) {
    double startx = x1+i * step;
    double endx = x1+(i+1)*step;

    double starty=f(p,startx);
    double endy = f(p,endx);

    double segmLength=Math.sqrt(step*step +(endy-starty)*(endy-starty));
    ans += segmLength;
}
         /////////////////// */
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
        double totalArea = 0;
        /** add you code below
         */
        ArrayList<Double> intersections = new ArrayList();
        intersections.add(x1);
        intersections.add(x2);
        int check = 100;
        double step = (x2 - x1) / check;
        for (int i = 0; i < check; i++) {
            double leftx = x1 + i * step;
            double rightx = x1 + (i + 1) * step;

            double difLeft = f(p1, leftx) - f(p2, leftx);
            double difRight = f(p1, rightx) - f(p2, rightx);

            if (difLeft * difRight <= 0) {
                double intersection = sameValue(p1, p2, leftx, rightx, EPS);


                boolean alrdeyExist = false;
                for (double exist : intersections) {
                    if (Math.abs(exist - intersection) <= EPS) {
                        alrdeyExist = true;
                        break;
                    }
                }
                if (!alrdeyExist) {
                    intersections.add(intersection);

                }
            }
        }
        Collections.sort(intersections);
        for (int i = 0; i < intersections.size() - 1; i++) {
            double subx1 = intersections.get(i);
            double subx2 = intersections.get(i + 1);
            double subArea = areaInsubRange(p1, p2, subx1, subx2, numberOfTrapezoid);
            totalArea += subArea;

        }
        return totalArea;
    }
private static  double areaInsubRange(double[]p1,double[]p2,double x1,double x2,int numberOfTrapezoid) {

        double area = 0;
        double width = (x2-x1)/numberOfTrapezoid;
for (int i=0;i<numberOfTrapezoid;i++) {
    double leftx = x1+i*width;
    double rightx = x1+(i+1)*width;

    double h1=Math.abs(f(p1,leftx)-f(p2,leftx));
    double h2=Math.abs(f(p1,rightx)-f(p2,rightx));
    area += width*(h1+h2)/2;
}
         /////////////////// */
		return area;
	}
	/**
	 * This function computes     the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
        double[] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below
         */
        if (p == null || p.trim().isEmpty() || p.equals("0")) return ZERO;
        p = p.replace(" ", "");
        int maxDegree = 0;
        String[] terms = p.split("(?=[-+])");
        for(String term : terms) {
            if(term.isEmpty()) continue;

            if(term.contains("x^")) {
                int caretPos = term.indexOf("^");
                int degree = Integer.parseInt(term.substring(caretPos + 1));
                maxDegree = Math.max(maxDegree, degree);
            } else if(term.contains("x")) {
                maxDegree = Math.max(maxDegree, 1);
            }
        }
        ans = new double[maxDegree + 1];

        for (String term : terms) {
                if (term.isEmpty()) continue;

                if (term.contains("x^")) {
                    // Term like "-1.0x^2" or "+3.5x^3"
                    int xPos = term.indexOf("x");
                    int caretPos = term.indexOf("^");

                    String coeffStr = term.substring(0, xPos);
                    if (coeffStr.equals("+") || coeffStr.equals("")) coeffStr = "1";
                    if (coeffStr.equals("-")) coeffStr = "-1";

                    double coeff = Double.parseDouble(coeffStr);
                    int degree = Integer.parseInt(term.substring(caretPos + 1));
                    ans[degree] = coeff;
                } else if (term.contains("x")) {
                    // Term like "3x" or "-2.0x"
                    int xPos = term.indexOf("x");
                    String coeffStr = term.substring(0, xPos);

                    if (coeffStr.equals("+") || coeffStr.equals("")) coeffStr = "1";
                    if (coeffStr.equals("-")) coeffStr = "-1";

                    double coeff = Double.parseDouble(coeffStr);
                    ans[1] = coeff;

                } else {
                    // Constant term like "5.0" or "-3"
                    double coeff = Double.parseDouble(term);
                    ans[0] = coeff;
                }
            }


        return ans;
        /////////////////// */
    }

	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
        /** add you code below
         */
int max=Math.max(p1.length,p2.length);
double[]ans = new double[max];
for( int i=0; i<p1.length;i++) {
    ans[i]=p1[i];
}
for( int i=0; i<p2.length;i++) {
    ans[i]+=p2[i];
}
         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below
         */
        int len=p1.length+p2.length-1;
        double[] mul = new double[len];
      if (p1.length==0||p2.length==0) return ans;
      for (int i=0;i<p1.length;i++) {
          for (int j=0;j<p2.length;j++) {
              mul[i+j]+=p1[i]*p2[j];
          }
      }
         /////////////////// */
		return mul;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
        /** add you code below
         */
if (po.length<=1) return ZERO;
double [] ans = new double[po.length-1];
for (int i=1;i<po.length;i++) {
    ans[i-1]=i*po[i];

}
         /////////////////// */
		return ans;
	}
}
