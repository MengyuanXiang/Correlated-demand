package multinormal.test;

import java.util.Arrays;


import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalGen;
import umontreal.iro.lecuyer.rng.MRG31k3p;
import umontreal.iro.lecuyer.rng.RandomStream;

public class Test {

	public static double[] multinormalPointGeneration(int N, double [] mu, double [][] sigma) {
		RandomStream s = new MRG31k3p();
		NormalGen standardNormal = new NormalGen(s, 0, 1);

		MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);

		double[] point = new double[N];
		gen.nextPoint(point);
		
		return point;
	
	}
	
	public static RealMatrix matrixInverse(RealMatrix m) {
		RealMatrix pInverse = new LUDecomposition(m).getSolver().getInverse();
		return pInverse;
		//System.out.println(pInverse.getEntry(0, 2));
		//System.out.println(pInverse.toString());
	}
	
	public static RealMatrix createReducedMatrix(RealMatrix m) {
		double[][] matrix = new double[m.getRowDimension()-1][m.getColumnDimension()-1];
		for(int i = 1; i < m.getRowDimension(); i++) {
			for(int j = 1; j < m.getColumnDimension(); j++) {
				matrix[i-1][j-1]=m.getEntry(i, j);
			}
		}
		RealMatrix n = MatrixUtils.createRealMatrix(matrix);
		return n;
		
	}



	public static RealMatrix matrix21(RealMatrix m) {
		return m.getSubMatrix(1, m.getRowDimension(), 0, 0);
	}
	
	

	public static RealMatrix computeConditionalExpectedDemand(double [] expDemand, double [] realizationDemand, RealMatrix m) {
		RealMatrix d1 =null;
		RealMatrix d2 =null;
		RealMatrix sigma21 =null;
		RealMatrix sigma11Inv =null;
		RealMatrix zeta1 =null;
		
		
		double[][] ed1 = new double[1][];
		ed1[0] = expDemand;
		RealMatrix ed1Matrix = MatrixUtils.createRealMatrix(ed1); 
		d1=ed1Matrix.getSubMatrix(0, 0, 0, 0).transpose(); //d1
		
	
		double[][] ed2 = new double[1][];
		ed2[0] = expDemand;
		RealMatrix ed2Matrix = MatrixUtils.createRealMatrix(ed2); 
		d2=ed2Matrix.getSubMatrix(0, 0, 1, expDemand.length - 1).transpose(); //d2
		
		sigma21 = m.getSubMatrix(1, m.getRowDimension() - 1, 0, 0);//sigma21
		
		sigma11Inv = matrixInverse(m.getSubMatrix(0, 0, 0, 0));//inverse of sigma11
		
		double[][] realisedDemand = new double[1][];
		realisedDemand[0] = realizationDemand;
		RealMatrix realisedDemandMatrix = MatrixUtils.createRealMatrix(realisedDemand); 
		zeta1 = realisedDemandMatrix.getSubMatrix(0, 0, 0, 0).transpose();//zeta1
		
		
		RealMatrix results = d2.add(sigma21.multiply(sigma11Inv.multiply(zeta1.subtract(d1))));		
		
        return results.transpose();	//conditional mean demand	
	}


	public static double[][] calucluteCovariance(double [] ExpDemand, double cv, double [] rho){
		double [] stdDemand =new double [ExpDemand.length];
		for (int i = 0; i < ExpDemand.length; i ++) {
			stdDemand[i] = cv*ExpDemand[i];
		}
		
		double [][] covariance = new double [ExpDemand.length][ExpDemand.length];
		
		//nunber of correlation is 1
		//for (int row=0; row<covariance.length;row++) {
		//	for (int col=0; col<covariance[row].length;col++) {
		//		if (row==col) {
		//			covariance[row][col]=stdDemand[row]*stdDemand[col];
		//		} else if (col==row+1 | col==row-1) {
		//			covariance[row][col]=stdDemand[row]*stdDemand[col]*rho;
		//		} else  {
		//			covariance[row][col]=0;
		//		}
		//	}
		//}
		
		//number of correlation is n
		for (int row=0; row<covariance.length;row++) {
			covariance[row][row]=stdDemand[row]*stdDemand[row];
			
			for (int nbcorrelation=0; nbcorrelation<rho.length; nbcorrelation++) {
				if (row+nbcorrelation<covariance.length-1) {
					covariance[row][row+nbcorrelation+1]=rho[nbcorrelation]*stdDemand[row]*stdDemand[row+nbcorrelation+1];
				}
				
				if (row-nbcorrelation>0) {
					covariance[row][row-nbcorrelation-1]=rho[nbcorrelation]*stdDemand[row]*stdDemand[row-nbcorrelation-1];
				}
			}
		}
		
		return covariance;//initial covariance matrix
		
	}

	public static void main(String[] args) {
		
		double[] ExpDemand = {5,8,24,39,16,29,51,39,75,69,26,20,32,11,19};
		double cv = 0.1;
		double[] rho = {0.5, 0.3, 0.1};
		
		double[][] covMatrix = calucluteCovariance(ExpDemand,cv,rho);
		RealMatrix covarianceMatrix = MatrixUtils.createRealMatrix(covMatrix); 
		System.out.println(covarianceMatrix.toString());
		
		double [] realDemand = multinormalPointGeneration(ExpDemand.length, ExpDemand, covMatrix);
		RealMatrix realDemandMatrix = MatrixUtils.createColumnRealMatrix(realDemand); 
		System.out.println(realDemandMatrix.toString());
		
		RealMatrix conditionalMeanDemand = computeConditionalExpectedDemand(ExpDemand,realDemand, covarianceMatrix); 
		System.out.println(conditionalMeanDemand.toString());
		
		RealMatrix M = MatrixUtils.createRealMatrix(covMatrix);
		
		RealMatrix inverseM =  matrixInverse(M);
		RealMatrix reducedInverseM = createReducedMatrix(inverseM);
		RealMatrix conditionalCovarianceMatrix =  matrixInverse(reducedInverseM);
		System.out.println(conditionalCovarianceMatrix.toString());

			
		 
	}

}
