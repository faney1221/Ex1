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
            Creates a polynomial from 2 points (line) or 3 points (parabola).
            If more or fewer points → returns null.
             */
            public static double[] PolynomFromPoints(double[] xx, double[] yy) {
                double [] ans =null;
                if ( xx==null || yy==null ) {
                    return null;
                }
                if (xx.length!=yy.length) { return null; }
                if (xx.length<2 || xx.length>3) { return null; }

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
                 double b =(y0-y1)/(x0-x1)-a*(x0+x1);
                 double c=y0-a*x0*x0-b*x0;
                 ans = new double[] {c,b,a};
                }
        
                return ans;
        
        
            }
            /** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
             * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
             * @param p1 first polynomial function
             * @param p2 second polynomial function
             * @return true iff p1 represents the same polynomial function as p2.

             * Checks if two polynomials represent the same function by comparing f(x) at several x values.
             */
            public static boolean equals(double[] p1, double[] p2) {
                boolean ans = true;
                /** add you code below
                 */

                    int maxDegree = Math.max(p1.length, p2.length);

                    for (int i = 0; i <= maxDegree; i++) {
                        double x = i;
                        if (Math.abs(f(p1, x) - f(p2, x)) > EPS) return false;
                    }
                    return true;
                }






                 /////////////////// */

        
            /** 
             * Computes a String representing the polynomial function.
             * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
             * @param poly the polynomial function represented as an array of doubles
             * @return String representing the polynomial function:
             *
             * Creates a string representation of the polynomial.
             * Example: {2, 0, 3} → "3x^2+2"
             */
            public static String poly(double[] poly) {
                String ans = "";
                if(poly.length==0) {ans="0";}
                else {
                    /** add you code below
                     */
        for (int i=poly.length-1;i>=0;i--){
            if(poly[i]==0.0||poly[i]==-0.0) continue;
            if ( !ans.isEmpty()&& poly[i]>0) ans+= "+";
        
            if (poly[i]<0) ans+= "";
        
            if (i==0){
                ans += poly[i];}
        
             else if(i==1){ans += poly[i] +"x";
            }
             else {ans += poly[i]+"x^"+i;}
            }
        }
                if (ans.isEmpty()) {
                    ans="0";
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
             *
             * Finds x where p1(x) ≈ p2(x).
             * It builds a difference polynomial p1 - p2 and calls root_rec.
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
             * @param n - (A positive integer value (1,2,...).
             * @return the length approximation of the function between f(x1) and f(x2).
             *
             * Approximates the arc length of the polynomial curve between x1 and x2.
             * Uses n straight-line segments:
             */
            public static double length(double[] p, double x1, double x2, int n) {
                /** add you code below
                 */

                    if (n <= 0) return 0;
                    if (x1 == x2) return 0;

                    if (x2 < x1) { double t = x1; x1 = x2; x2 = t; }

                    double dx = (x2 - x1) / n;
                    double sum = 0;

                    double prevX = x1;
                    double prevY = f(p, prevX);

                    for (int i = 1; i <= n; i++) {
                        double currX = x1 + i * dx;
                        double currY = f(p, currX);

                        double dxSeg = currX - prevX;
                        double dySeg = currY - prevY;

                        sum += Math.sqrt(dxSeg*dxSeg + dySeg*dySeg);

                        prevX = currX;
                        prevY = currY;
                    }

                    return sum;
                }





                /////////////////// */

            
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
             *
             * Computes the area between two polynomials using trapezoids.
             * Steps:
             * Find intersection points.
             * Split the range.
             * Compute area in each sub-range using trapezoids
             */
            public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
                double totalArea = 0;
                /** add you code below
                 */
                    if (x2 < x1) { double t = x1; x1 = x2; x2 = t; }

                    ArrayList<Double> xs = new ArrayList<>();
                    xs.add(x1);
                    xs.add(x2);

                    // Build diff polynomial
                    int m = Math.max(p1.length, p2.length);
                    double[] diff = new double[m];
                    for (int i = 0; i < m; i++) {
                        double a = i < p1.length ? p1[i] : 0;
                        double b = i < p2.length ? p2[i] : 0;
                        diff[i] = a - b;
                    }

                    // Detect intersections
                    int samples = 1000;
                    double step = (x2 - x1) / samples;
                    double prevX = x1;
                    double prevF = f(diff, prevX);

                    for (int i = 1; i <= samples; i++) {
                        double currX = x1 + i * step;
                        double currF = f(diff, currX);

                        if (prevF * currF <= 0) {
                            double root = sameValue(p1, p2, prevX, currX, EPS);

                            boolean exists = false;
                            for (double v : xs) {
                                if (Math.abs(v - root) < EPS) { exists = true; break; }
                            }
                            if (!exists) xs.add(root);
                        }

                        prevX = currX;
                        prevF = currF;
                    }

                    // Sort intervals
                    Collections.sort(xs);

                    // Compute area in each interval
                    double total = 0;
                    for (int i = 0; i < xs.size() - 1; i++) {
                        double L = xs.get(i);
                        double R = xs.get(i + 1);
                        total += areaInsubRange(p1, p2, L, R, numberOfTrapezoid);
                    }

                    return total;
                }


        public static  double areaInsubRange(double[]p1,double[]p2,double x1,double x2,int numberOfTrapezoid) {

                double dx = (x2 - x1) / numberOfTrapezoid;
                double sum = 0;

                for (int i = 0; i < numberOfTrapezoid; i++) {
                    double lx = x1 + i * dx;
                    double rx = lx + dx;
                    double h1 = Math.abs(f(p1, lx) - f(p2, lx));
                    double h2 = Math.abs(f(p1, rx) - f(p2, rx));
                    sum += dx * (h1 + h2) / 2;
                }
                return sum;
            }


                 /////////////////// */


            /**
             * This function computes     the array representation of a polynomial function from a String
             * representation. Note:given a polynomial function represented as a double array,
             * getPolynomFromString(poly(p)) should return an array equals to p.
             * 
             * @param p - a String representing polynomial function.
             *
             *          arses a polynomial string like:
             * "-1.2x^3 + 3.1x^2 + 2"
             * and returns the coefficient array.
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
             *
             * Adds two polynomials term by term.
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
             *
             * Multiplies two polynomials using nested loops.
             * @return
             */
            public static double[] mul(double[] p1, double[] p2) {
                /** add you code below
                 */
                if (p1.length == 0||p2.length==0) return ZERO;
                int len=p1.length+p2.length-1;
                double[] mul = new double[len];
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
             * eturns the derivative polynomial.
             * Example:
             * {2, 3, 4} → derivative = {3, 8} (3 + 8x)
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