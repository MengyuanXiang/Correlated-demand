package simulation;



import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Scanner;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.mahout.math.Arrays;

import sample_opl.Model;
import umontreal.iro.lecuyer.randvar.NormalGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalCholeskyGen;
import umontreal.iro.lecuyer.randvarmulti.MultinormalGen;
import umontreal.iro.lecuyer.rng.MRG31k3p;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.stat.Tally;

public class Simulator {
	
	public static double[][] calculateCovariance(double [] ExpDemand, double cv, double rho){
		double [] stdDemand =new double [ExpDemand.length];
		for (int i = 0; i < ExpDemand.length; i ++) {
			stdDemand[i] = cv*ExpDemand[i];
		}
		
		double [][] covariance = new double [ExpDemand.length][ExpDemand.length];
		
		
		for (int row=0; row<covariance.length;row++) {
			for (int col=0; col<covariance[row].length;col++) {
				if (row==col) {
					covariance[row][col]=stdDemand[row]*stdDemand[col];
				} else if (col==row+1 | col==row-1) {
					covariance[row][col]=stdDemand[row]*stdDemand[col]*rho;
				} else  {
					covariance[row][col]=0;
				}
			}
		}
		
		/*
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
		}*/
		
		return covariance;		
	}
	
	public static RealMatrix matrixInverse(RealMatrix m) {
		RealMatrix pInverse = new LUDecomposition(m).getSolver().getInverse();
		return pInverse;
		//System.out.println(pInverse.getEntry(0, 2));
		//System.out.println(pInverse.toString());
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
		
		sigma21 = m.getSubMatrix(1, m.getRowDimension() - 1, 0, 0);
		
		sigma11Inv = matrixInverse(m.getSubMatrix(0, 0, 0, 0));
		
		double[][] realisedDemand = new double[1][];
		realisedDemand[0] = realizationDemand;
		RealMatrix realisedDemandMatrix = MatrixUtils.createRealMatrix(realisedDemand); 
		zeta1 = realisedDemandMatrix.getSubMatrix(0, 0, 0, 0).transpose();
		
		
		RealMatrix results = d2.add(sigma21.multiply(sigma11Inv.multiply(zeta1.subtract(d1))));		
		
        return results.transpose();		
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
	
	public static double[] simulateOnePeriod(int t, 
										   double[] realizations,
										   double initialStock,
										   double ordercost,
										   double holdingcost,
										   double penaltycost,
										   double unitcost,
										   double[] expDemand,
										   double cv,
										   double rho,
										   int Nbpartitions){
		double[] previousRealizations = new double[t];
		System.arraycopy(realizations, 0, previousRealizations, 0, t);
		
		double[][] covariance = calculateCovariance(expDemand, cv, rho);
		
		int Nbperiods = expDemand.length-t;
		
		double[] shortexpDemand = expDemand;
		RealMatrix conditionalCovarianceMatrix = MatrixUtils.createRealMatrix(covariance);
		RealMatrix M = conditionalCovarianceMatrix;
		double[] realizationsBuffer = realizations;
		for(int k = 0; k < t; k++) {
			RealMatrix inverseM =  matrixInverse(conditionalCovarianceMatrix);
			RealMatrix reducedInverseM = createReducedMatrix(inverseM);
			conditionalCovarianceMatrix =  matrixInverse(reducedInverseM);
			shortexpDemand = computeConditionalExpectedDemand(shortexpDemand, realizationsBuffer, M).getRow(0);
			double[][] realizationMatrix = new double[1][];
			realizationMatrix[0] = realizationsBuffer;
			RealMatrix reducedRealizations = MatrixUtils.createRealMatrix(realizationMatrix);
			realizationsBuffer =  reducedRealizations.getSubMatrix(0, 0, 1, realizationsBuffer.length - 1).getRow(0);
			M = conditionalCovarianceMatrix;
		}
		double [][] shortCovariance = conditionalCovarianceMatrix.getData();
		
		double Q = Double.NaN;
		try{
	         Model model = new Model(
	        	   Nbperiods, 
	        	   shortexpDemand, 
	               ordercost, 
	               holdingcost,
	               penaltycost,
	               unitcost, 
	               initialStock, 
	               shortCovariance,
	               Nbpartitions,
	               null);
	         Q = model.solve("rs_milp_piecewise_penalty_ub") - initialStock;
	         
	         double epsilon = 0.1;
	         if(Q < epsilon) Q = 0;
	         
	      }catch(Exception e){
	         e.printStackTrace();
	      }
		
		
		double hc = holdingcost*Math.max(initialStock + Q - realizations[t], 0);
		double pc = penaltycost*Math.max(realizations[t] - (initialStock + Q), 0);
		double oc = unitcost*Q + (Q>0?ordercost:0);
		return new double[] {hc+pc+oc, initialStock + Q - realizations[t]};
	}
	
	static RandomStream rng = new MRG31k3p();
	
	public static double simulateOneRun(
			   double[] realizations,
			   double initialStock,
			   double ordercost,
			   double holdingcost,
			   double penaltycost,
			   double unitcost,
			   double[] expDemand,
			   double cv,
			   double rho,
			   int Nbpartitions) {
		double etc = 0;
		for(int t = 0; t < realizations.length; t++) {
			double[] result = simulateOnePeriod(t,
					                            realizations, 
					                            initialStock, 
					                            ordercost, 
					                            holdingcost, 
					                            penaltycost, 
					                            unitcost,
					                            expDemand,
					                            cv,
					                            rho,
					                            Nbpartitions);
			initialStock = result[1];
			etc += result[0];
			//System.out.println(result[0] + "\t" + result[1]);
		}
		//System.out.println(etc);
		return etc;
	}

	/*
	public static void main(String[] args) throws FileNotFoundException {
		double penaltycost = 10 ;
		double unitcost = 0;
		double ordercost = 100;
		double[] expDemand = {5,8,24,39,16,29,51,39,75,69,26,20,32,11,19};
		double holdingcost = 1;
		double initialStock = 0;
		int Nbpartitions = 10;
		
		double cv = 0.1;
		double rho = 0.5;
		
		
		String fileName="/Users/mengyuanxiang/Downloads/myFile32.txt";
		
		File textFile = new File(fileName);
		
		Scanner in = new Scanner(textFile);
		
		int rows = 5000;//number of simulations
		int columns = 15;//number of time perios
		double[][] realizations = new double [rows][columns];
		
		while(in.hasNextLine()) {
			for (int i=0; i<realizations.length;i++) {
				String[] line = in.nextLine().trim().split(" ");
				for (int j=0;j<line.length;j++) {
					realizations[i][j]=Double.parseDouble(line[j]);
					
				}
			}
			//String line = in.nextLine();
			//System.out.println(line);
			
		}
		//System.out.println(Arrays.toString(realizations));
		in.close();*/

		/*
		double[][] realizations = {
				{14.8390050967702,19.1946573050903,10.8685800582842,14.1260021324616,16.3334884846290,28.0399277361043,17.4013636163451,51.5622851073546},
				{10.6470943209089,20.9752720362035,12.6841375336190,11.4111973550090,17.0291562340996,30.6709684851637,16.1195622474581,52.4440377413956},
				{13.3790442329951,19.6703767978565,10.1858678943706,10.1613882338950,19.4695962446614,37.5198618262082,18.1243971842410,46.1534503026573},
				{10.3017615885902,15.0023167265418,10.7992738466370,10.8204497563823,20.8604147177837,28.0641911222693,16.0897395844909,59.0349021312539},
				{10.5488915698347,13.1496738060104,9.36289107505558,10.8381733603603,17.4128902947308,37.7552301114486,20.8337091381797,55.8344712550846}	
		};
		*/
/*
		double[] results = new double[realizations.length];
		for(int i = 0; i< realizations.length; i++) {
			double etc = simulateOneRun(realizations[i], initialStock, ordercost, holdingcost,
				penaltycost, unitcost, expDemand, cv, rho, Nbpartitions);
			results[i] = etc;
		}
		
		double totalcost=0;
		for(int j =0; j<results.length;j++) {
			totalcost += results[j];
		}
		//System.out.println(Arrays.toString(results));
		System.out.println(totalcost);
		System.out.println(totalcost/results.length);
	}
	*/
	
	
	public static double[] multinormalPointGeneration(int N, double [] mu, double [][] sigma) {
		
		NormalGen standardNormal = new NormalGen(rng, 0, 1);
		MultinormalGen gen = new MultinormalCholeskyGen(standardNormal, mu, sigma);

		double[] point = new double[N];
		gen.nextPoint(point);
		return point;	
		
	}

	public static void main(String[] args) {
		int counter=0;
		File file = new File("results.txt");
		
		double errorThreshold = 0.005;
		double confidenceProbability = 0.95;

		
		try {
			java.text.NumberFormat format = new java.text.DecimalFormat("0.01");
			PrintWriter pw = new PrintWriter(new FileWriter(file));

			double[] penaltyCostArray = {5,10};
			double[] orderCostArray = { 200, 400 };
			double[] rhoArray = { -0.5,-0.25,0.25, 0.5 };
			double[] cvArray = {0.1,0.2};

			for (int rh = 0; rh < rhoArray.length; rh++) {
				double rho = rhoArray[rh];

				for (int o = 0; o < orderCostArray.length; o++) {
					double ordercost = orderCostArray[o];

						for (int p = 0; p < penaltyCostArray.length; p++) {
                      	    double penaltycost = penaltyCostArray[p];
                      	    
                      	    for (int c = 0; c< cvArray.length; c++) {
                      	    	double cv = cvArray[c];
                      	  


							double[] expDemand = {105,46,17,40,22,38,55,58,112,11,24,68,7,104,43}; 
							double unitcost =0;
							double holdingcost = 1;
							double initialStock = 0;
							int Nbpartitions = 10;

							Tally tally = new Tally();
							rng.resetStartStream(); 
							boolean stop = false;
							for (int r = 0; r < 20 || !stop; r++) {
								// rng.resetNextSubstream();
								
								
								double[] realizations = multinormalPointGeneration(expDemand.length, expDemand,
										calculateCovariance(expDemand, cv, rho));
										
								
								// double[] realizations = expDemand;
								double etc = simulateOneRun(realizations, initialStock, ordercost, holdingcost,
										penaltycost, unitcost, expDemand, cv, rho, Nbpartitions);
								
								tally.add(etc);
								
								System.out.print(".");
								if(r >= 20) {
									double[] centerAndRadius = new double[2];
									tally.confidenceIntervalStudent(confidenceProbability, centerAndRadius);
									if(r%10==0) {
										System.out.println(tally.report());
										System.out.println("Current error: "+centerAndRadius[1]/centerAndRadius[0]);
									}
									if(centerAndRadius[1]/centerAndRadius[0] <= errorThreshold) {
										stop = true;
										System.out.println(r);
									}
								}
							}

							counter++;
							System.out.println("********************************************");
							System.out.println(counter);
							System.out.println(tally.report());
							System.out.println("********************************************");
							pw.println(format.format(tally.average()));
							
						}
						
					}
				
				}
			}
			pw.close();

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
