/*
 * PROJECT III: TriMatrix.java
 *
 * This file contains a template for the class TriMatrix. Not all methods are
 * implemented. Make sure you have carefully read the project formulation
 * before starting to work on this file. You will also need to have completed
 * the Matrix class.
 *
 * Remember not to change the names, parameters or return types of any
 * variables in this file!
 *
 * The function of the methods and instance variables are outlined in the
 * comments directly above them.
 */

public class TriMatrix extends Matrix {
    /**
     * An array holding the diagonal elements of the matrix.
     */
    private double[] diag;

    /**
     * An array holding the upper-diagonal elements of the matrix.
     */
    private double[] upper;

    /**
     * An array holding the lower-diagonal elements of the matrix.
     */
    private double[] lower;
    
    /**
     * Constructor function: should initialise m and n through the Matrix
     * constructor and set up the data array.
     *
     * @param N  The dimension of the array.
     */
    public TriMatrix(int N) {
        // You need to fill in this method.
    	super(N,N);
    	if (N <= 0) {
    		throw new MatrixException("Invalid matrix dimensions.");
    	}
    	this.diag = new double[N];
    	this.upper = new double[N-1];
    	this.lower = new double[N-1];
    }
    
    /**
     * Getter function: return the (i,j)'th entry of the matrix.
     *
     * @param i  The location in the first co-ordinate.
     * @param j  The location in the second co-ordinate.
     * @return   The (i,j)'th entry of the matrix.
     */
    public double getIJ(int i, int j) {
        // You need to fill in this method.
    	if (i < 0 || i >= this.m) {
    		throw new MatrixException("The matrix index position is out of bounds");
    	}
    	else if (j < 0 || j >= this.n) {
    		throw new MatrixException("The matrix index position is out of bounds");
    	}
    	
    	if (Math.abs(i-j) > 1) {
    		return (double) 0;
    	}
    	else {
    		if (i == j) {
    			return this.diag[i];
    		}
    		else if (i > j) {
    			return this.lower[j];
    		}
    		else {
    			return this.upper[i];
    		}
    	}
    }
    
    /**
     * Setter function: set the (i,j)'th entry of the data array.
     *
     * @param i    The location in the first co-ordinate.
     * @param j    The location in the second co-ordinate.
     * @param val  The value to set the (i,j)'th entry to.
     */
    public void setIJ(int i, int j, double val) {
        // You need to fill in this method.
    	if (i < 0 || i >= this.m) {
    		throw new MatrixException("The matrix index position is out of bounds");
    	}
    	else if (j < 0 || j >= this.n) {
    		throw new MatrixException("The matrix index position is out of bounds");
    	}
    	else if (Math.abs(i - j) > 1) {
    		throw new MatrixException("The tri-matrix index position is out of bounds");
    	}
    	
    	if (i == j) {
			this.diag[i] = val;
		}
		else if (i > j) {
			this.lower[j] = val;
		}
		else if (j > i){
			this.upper[i] = val;
		}
	
    }
    
    /**
     * Return the determinant of this matrix.
     *
     * @return The determinant of the matrix.
     */
    public double determinant() {
        // You need to fill in this method.
    	TriMatrix dMatrix = this.decomp();
    	double det = dMatrix.getIJ(0, 0);
    	
    	
    	for (int i = 1; i < this.n; i++) {
    		det = det*dMatrix.getIJ(i, i);
    	}
    	
    	return det;
    }
    
    /**
     * Returns the LU decomposition of this matrix. See the formulation for a
     * more detailed description.
     * 
     * @return The LU decomposition of this matrix.
     */
    public TriMatrix decomp() {
        // You need to fill in this method.
    	TriMatrix dMatrix = new TriMatrix(this.n);
    	double dstar, lstar, ustar;
    	
    	if (this.getIJ(0, 0) == 0) {
    		throw new MatrixException("The matrix is singular");
    	}
    	dstar = this.getIJ(0, 0);
    	dMatrix.setIJ(0, 0, dstar);
    	dMatrix.setIJ(0, 1, this.getIJ(0, 1));
    	for (int i = 2; i < this.n; i++) {
    		if (this.getIJ(0, 0) == 0) {
        		throw new MatrixException("The matrix is singular");
        	}
    		ustar = this.getIJ(i-1, i);
    		lstar = this.getIJ(i, i-1)/dstar;
    		dstar = this.getIJ(i-1, i-1) - ((this.getIJ(i-1, i-2)*this.getIJ(i-2, i-1))/dstar);
    		
    		dMatrix.setIJ(i-1, i, ustar);
    		dMatrix.setIJ(i-1, i-2, lstar);
    		dMatrix.setIJ(i-1, i-1, dstar);
    		
    	}
    	if (this.getIJ(0, 0) == 0) {
    		throw new MatrixException("The matrix is singular");
    	}
		lstar = this.getIJ(this.n-1, this.n-2)/dstar;
		dstar = this.getIJ(this.n-1, this.n-1) - ((this.getIJ(this.n-1, this.n-2)*this.getIJ(this.n-2, this.n-1))/dstar);
		
		dMatrix.setIJ(this.n-1, this.n-2, lstar);
		dMatrix.setIJ(this.n-1, this.n-1, dstar);
		
		return dMatrix;
    	
    }

    /**
     * Add the matrix to another matrix A.
     *
     * @param A  The Matrix to add to this matrix.
     * @return   The sum of this matrix with the matrix A.
     */
    public Matrix add(Matrix A){
        // You need to fill in this method.
    	if (this.n != A.m || this.n != A.n) {
    		throw new MatrixException("Invalid matrix dimensions for addition");
    	}
    	
    	Matrix sumA = new GeneralMatrix((GeneralMatrix) A);
    	for (int i = 0; i < this.n - 1; i++) {
    		sumA.setIJ(i, i, this.diag[i] + A.getIJ(i, i));
    		sumA.setIJ(i, i + 1, this.upper[i] + A.getIJ(i, i + 1));
    		sumA.setIJ(i + 1, i, this.lower[i] + A.getIJ(i + 1, i));
    	}
    	sumA.setIJ(this.n - 1, this.n - 1, this.diag[this.n - 1] + A.getIJ(this.n - 1, this.n - 1));
    	
    	return sumA;
    }
    
    /**
     * Multiply the matrix by another matrix A. This is a _left_ product,
     * i.e. if this matrix is called B then it calculates the product BA.
     *
     * @param A  The Matrix to multiply by.
     * @return   The product of this matrix with the matrix A.
     */
    public Matrix multiply(Matrix A) {
        // You need to fill in this method.
    	if (this.n != A.m) {
    		throw new MatrixException("Incorrect matrix dimensions for multiplication");
    	}
    	
    	Matrix multMatrix = new GeneralMatrix(this.n, A.n);
    	double term;
    	
    	for (int j = 0; j < A.n; j++) {
    		term = this.diag[0]*A.getIJ(0, j) + this.upper[0]*A.getIJ(1, j);
    		multMatrix.setIJ(0, j, term);
    	}
    	
    	for (int i = 1; i < this.n-1; i++) {
    		for (int j = 0; j < A.n; j++) {
    			term = this.lower[i-1]*A.getIJ(i-1, j) + this.diag[i]*A.getIJ(i, j) + this.upper[i]*A.getIJ(i+1, j);
    			multMatrix.setIJ(i, j, term);
    		}
    	}
    	
    	for (int j = 0; j < A.n; j++) {
    		term = this.lower[this.n-2]*A.getIJ(A.m-2, j) + this.diag[this.n-1]*A.getIJ(A.m-1, j);
    		multMatrix.setIJ(this.n-1, j, term);
    	}
    	
    	return multMatrix;
    	
    	
    	
    }
    
    /**
     * Multiply the matrix by a scalar.
     *
     * @param a  The scalar to multiply the matrix by.
     * @return   The product of this matrix with the scalar a.
     */
    public Matrix multiply(double a) {
        // You need to fill in this method.
    	Matrix scalMatrix = new GeneralMatrix(this.n, this.n);
    	for (int i = 0; i < this.n - 1; i++) {
    		scalMatrix.setIJ(i, i, this.diag[i]*a);
    		scalMatrix.setIJ(i, i+1, this.upper[i]*a);
    		scalMatrix.setIJ(i+1, i, this.lower[i]*a);
    	}
    	scalMatrix.setIJ(this.n-1, this.n-1, this.diag[this.n - 1]*a);
    	
    	return scalMatrix;
    }

    /**
     * Populates the matrix with random numbers which are uniformly
     * distributed between 0 and 1.
     */
    public void random() {
        // You need to fill in this method.
    	for (int i = 0; i < this.n-1; i++) {
    		this.diag[i] = Math.random();
    		this.lower[i] = Math.random();
    		this.upper[i] = Math.random();
    	}
    	this.diag[this.n-1] = Math.random();
    }
    
    /*
     * Your tester function should go here.
     */
    public static void main(String[] args) {
        // You need to fill in this method.
    	TriMatrix A = new TriMatrix(3);
    	GeneralMatrix B = new GeneralMatrix(3,3);
    	
    	B.random();
    	A.random();
    	
    	System.out.println("A = \n" + A.toString() + "\n");
    	System.out.println("B = \n" + B.toString() + "\n");
    	
    	System.out.println("AB = \n" + A.multiply(B) + "\n");
    	System.out.println("9*A = \n" + A.multiply(9) + "\n");
    	
    	System.out.println("A + B = \n" + A.add(B) + "\n");
    	System.out.println("det(A) = \n" + A.determinant() + "\n");
    }
}