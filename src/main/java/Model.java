

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;

import org.apache.commons.math3.distribution.PoissonDistribution;

import ilog.concert.IloException;
import ilog.opl.IloCplex;
import ilog.opl.IloCustomOplDataSource;
import ilog.opl.IloOplDataHandler;
import ilog.opl.IloOplDataSource;
import ilog.opl.IloOplErrorHandler;
import ilog.opl.IloOplFactory;
import ilog.opl.IloOplModel;
import ilog.opl.IloOplModelDefinition;
import ilog.opl.IloOplModelSource;
import ilog.opl.IloOplSettings;

public class Model {
	int nbitems;
	int nbmonths;
	double[] expDemand; 
	double groupfc;
	double[] ordercost;
	double holdingcost; 
	double penaltycost;
	int[] leadtime;
	int[] initialStock;
	double[][] cyclecost;
	int[][] P;
	
	String instanceIdentifier;


	public Model(
			int nbitems,
			int nbmonths, 
			double groupfc,
			double[][] cyclecost,
			int[][] P,
			String instanceIdentifier
			){
		this.nbitems = nbitems;
		this.nbmonths = nbmonths;
		this.groupfc = groupfc;
		this.cyclecost =cyclecost;
		this.P = P;
		this.instanceIdentifier = instanceIdentifier;
	}
	
	   private InputStream getMILPModelStream(File file){
		      FileInputStream is = null;
		      try{
		         is = new FileInputStream(file);
		      }catch(IOException e){
		         e.printStackTrace();
		      }
		      return is;
		   }
	
	// calculate order-up-to-position S for all items
	public static int[][] calculateS(
			int nbitems, 
			int nbmonths, 
			double holdingcost, 
			double penaltycost, 
			double [] expDemand,  
			int [] leadtime){
		int[][] S = new int[nbitems][nbmonths];
		for(int n = 0; n<nbitems;n++) {
			double [] lambda = new double[nbmonths];
			for (int t =0; t<nbmonths; t++) {			
				lambda[t]=(t+leadtime[n]+1)*expDemand[n];	
				//System.out.print(lambda[t]+"   ");
			}

			for (int t =0; t<nbmonths; t++) {
				for (int r = 1; r <10000; r++) {
					double lefthandside = 0;
					for (int i =0; i<t+1; i++) {
						PoissonDistribution pd = new PoissonDistribution(lambda[i]);
						lefthandside = pd.cumulativeProbability(r)+lefthandside;
					}

					if (lefthandside>=(t+1)*penaltycost/(penaltycost+holdingcost)) {
						S[n][t] = r;
						break;
					}
				}
			}
		}
		//System.out.println(S.toString());
		return S;
	}
	
	
	//calculate expected holding cost
	public static double calculateExpHoldingcost(double holdingcost, double meandemand, int x) {
		double expHolding = 0;
		PoissonDistribution pd = new PoissonDistribution(meandemand);
		for (int d=0; d<10000;d++) {
			expHolding =holdingcost*Math.max(x-d, 0)*pd.probability(d)+expHolding;
		}
		return expHolding;
	}

	//calculate expected penalty cost
	public static double calculateExpPenaltycost(double penaltycost, double meandemand, int x) {
		double exppenalty = 0;
		PoissonDistribution pd = new PoissonDistribution(meandemand);
		for (int d=0; d<10000;d++) {
			exppenalty =penaltycost*Math.max(d-x, 0)*pd.probability(d)+exppenalty;
		}
		return exppenalty;
	}
	
	//calculate the first period cyclecost without issuing an order	
	public static double[][] calculateFirstPeriodcostwithout(
			int nbitems, 
			int nbmonths, 
			double holdingcost, 
			double penaltycost, 
			double [] expDemand, 
			int [] initialStock){
		double[][] costWithout =new double[nbitems][nbmonths];
		for(int n =0; n<costWithout.length; n++) {
			for(int j =0;j<costWithout[n].length;j++) {
				double immediatecostWithout = 0;
				for(int i=0; i<j+1;i++) {
					double hwithout = calculateExpHoldingcost(holdingcost,(i+1)*expDemand[n],initialStock[n]);
					double pwithout = calculateExpPenaltycost(penaltycost,(i+1)*expDemand[n],initialStock[n]);
					immediatecostWithout = hwithout+pwithout+immediatecostWithout;
				}
				costWithout[n][j]=immediatecostWithout;
			}
		}
		return costWithout;
	}
	
	
	//calculate the first period cyclecost with issuing an order
	public static double[][] calculateFirstPeriodcostwithOrder(
			int [][] S, 
			int nbitems, 
			int nbmonths, 
			double holdingcost, 
			double penaltycost, 
			double [] ordercost, 
			double [] expDemand, 
			int [] leadtime){
		double[][] costWithOrder=new double[nbitems][nbmonths*nbmonths];
		for(int n =0; n<nbitems; n++) {
			int [] trueS = new int[nbmonths];
			for (int t =0; t<nbmonths; t++) {
				trueS[t] = (int)(S[n][t]-leadtime[n]*expDemand[n]);
			}
			System.out.print(Arrays.toString(trueS));
			
			for (int i = 0; i<nbmonths;i++) {
				for (int j = i;j<nbmonths;j++) {
					double immediatecostWith = 0;
					for (int m =1; m<j-i+1+1;m++) {
						double hwith = calculateExpHoldingcost(holdingcost,m*expDemand[n],trueS[j-i]);
					    double pwith = calculateExpPenaltycost(penaltycost,m*expDemand[n],trueS[j-i]);
						immediatecostWith = hwith+pwith+immediatecostWith;
					}
					costWithOrder[n][i*nbmonths+j] = immediatecostWith+ordercost[n];
				}
			}
		}
			
		return costWithOrder;
	}

	//calculate cyclecost
	public static double[][] calculateCyclecost(int nbitems, int nbmonths, double[][] costWithout, double[][] costWithOrder){
		double[][] cyclecost = new double[nbitems][nbmonths*nbmonths];
		cyclecost = costWithOrder;
		for (int n =0;n<nbitems; n++) {
			for (int t=0; t<nbmonths;t++) {
				if(costWithout[n][t]<=costWithOrder[n][t]) {
					cyclecost[n][t] = costWithout[n][t];
				}else {
					cyclecost[n][t] = costWithOrder[n][t];
				}	
			}
		}
		return cyclecost;
	}
	
	//record first period order decision
	public static int[][] firstPeriodOrderDecisin(int nbitems, int nbmonths, double[][] costWithout, double[][] costWithOrder){
		int[][] firstPeriodDecision = new int[nbitems][nbmonths];
		for (int n =0;n<nbitems; n++) {
			for (int t=0; t<nbmonths;t++) {
				if(costWithout[n][t]<=costWithOrder[n][t]) {
					firstPeriodDecision[n][t] = 0;
				}else {
					firstPeriodDecision[n][t] = 1;
				}	
			}
		}
		return firstPeriodDecision;
	}
	

	public static void main(String[] args) {
		int nbitems = 2;
        int nbmonths = 5; 
        double[] expDemand = {2, 3}; 
        double groupfc = 200;
        double[] ordercost = {80, 80}; 
        double holdingcost = 100; 
        double penaltycost = 2000;
        int[] leadtime = {0, 1};
        int[] initialStock = {0, 0};
        
        int[][] S = calculateS(nbitems,nbmonths,holdingcost,penaltycost,expDemand,leadtime);
        System.out.println(Arrays.deepToString(S));
        
        double[][] costWithout = calculateFirstPeriodcostwithout(nbitems, nbmonths, holdingcost, penaltycost, expDemand, initialStock);
        System.out.println(Arrays.deepToString(costWithout));
        
        double[][] costwithorder = calculateFirstPeriodcostwithOrder(S, nbitems, nbmonths, holdingcost, penaltycost, ordercost, expDemand, leadtime);
        System.out.println(Arrays.deepToString(costwithorder));
  
        double[][] cyclecost= calculateCyclecost(nbitems,nbmonths,costWithout,costwithorder);
        System.out.println(Arrays.deepToString(cyclecost));
        
        int[][] P = firstPeriodOrderDecisin(nbitems,nbmonths,costWithout,costwithorder);
        System.out.println(Arrays.deepToString(P));
        
        double shortestPath = Double.NaN;
        try{
        	Model model = new Model(
        			nbmonths, 
        			nbitems,
        			groupfc,
        			cyclecost,
        			P,
        			null
        			);
            shortestPath = model.solve("shortestpath");
            
         }catch(IloException e){
            e.printStackTrace();
         }
	}
        
        public double solve(String model_name) throws IloException{
            //IloOplFactory.setDebugMode(true);
            IloOplFactory oplF = new IloOplFactory();
            IloOplErrorHandler errHandler = oplF.createOplErrorHandler(System.out);
            IloCplex cplex = oplF.createCplex();
            IloOplModelSource modelSource=oplF.createOplModelSourceFromStream(getMILPModelStream(new File("./opl_models/shortestpath/"+model_name+".mod")),model_name);
            IloOplSettings settings = oplF.createOplSettings(errHandler);
            IloOplModelDefinition def=oplF.createOplModelDefinition(modelSource,settings);
            IloOplModel opl=oplF.createOplModel(def,cplex);
            cplex.setParam(IloCplex.IntParam.Threads, 8);
            cplex.setParam(IloCplex.IntParam.MIPDisplay, 2);
            /*cplex.setParam(IloCplex.IntParam.VarSel, 1);
            cplex.setParam(IloCplex.IntParam.ZeroHalfCuts, 2);
            cplex.setParam(IloCplex.IntParam.ImplBd, 2);
            cplex.setParam(IloCplex.IntParam.FracCuts, 2);
            cplex.setParam(IloCplex.IntParam.GUBCovers, 2);
            cplex.setParam(IloCplex.IntParam.DisjCuts, 2);
            cplex.setParam(IloCplex.IntParam.Covers, 2);
            cplex.setParam(IloCplex.IntParam.Cliques, 2);
            cplex.setParam(IloCplex.IntParam.FlowCovers, 2);
            cplex.setParam(IloCplex.IntParam.FlowPaths, 2);
            cplex.setParam(IloCplex.IntParam.MIRCuts, 2);
            cplex.setParam(IloCplex.IntParam.MIPEmphasis, 3);
            */

            IloOplDataSource dataSource = new Model.MyData(oplF);
            opl.addDataSource(dataSource);
            opl.generate();

            cplex.setOut(null);
            
            double start = cplex.getCplexImpl().getCplexTime();
            boolean status =  cplex.solve();
            double end = cplex.getCplexImpl().getCplexTime();
            if ( status )
            {   
             double objective = cplex.getObjValue();
             double time = end - start;
                //System.out.println("OBJECTIVE: " + objective);  
                //s = new double[Nbmonths];
                int[] groupPurchase = new int[nbmonths];
                //boolean[][] itemPurchase = new boolean[nbitems][nbmonths];
                for(int i = 0; i < nbmonths; i++){
                   //s[i] = cplex.getValue(opl.getElement("sValue").asNumVarMap().get(1+i));
                   groupPurchase[i] = (int) Math.round(cplex.getValue(opl.getElement("grouppurchase").asIntVarMap().get(1+i)));
                   //S[i] = cplex.getValue(opl.getElement("stock").asNumVarMap().get(1+i))+expDemand[i];
                   //System.out.println("S["+(i+1)+"]="+S[i]);
                }
                opl.postProcess();
                //opl.printSolution(System.out);
                //opl.end();
                oplF.end();
                //errHandler.end();
                //cplex.end();
                System.gc();

                //return objective;
                //System.out.println(S[0]);
                return groupPurchase[0];
            } else {
                System.out.println("No solution!");
                //opl.end();
                oplF.end();
                //errHandler.end();
                //cplex.end();
                System.gc();
                return Double.NaN;
            } 
            
        }
       
       
        class MyData extends IloCustomOplDataSource
        {
        	MyData(IloOplFactory oplF)
        	{
        		super(oplF);
        	}

        	public void customRead()
        	{
        		IloOplDataHandler handler = getDataHandler();

        		handler.startElement("nbmonths");
        		handler.addIntItem(nbmonths);
        		handler.endElement();

        		handler.startElement("nbitems");
        		handler.addIntItem(nbitems);
        		handler.endElement();

        		handler.startElement("groupfc");
        		handler.addNumItem(groupfc);
        		handler.endElement();

        		handler.startElement("cyclecost");
        		handler.startArray();
        		for (int i = 0 ; i<cyclecost.length; i++) {
        			handler.startArray();
        			for (int j = 0 ; j<cyclecost[i].length; j++) 
        				handler.addNumItem(cyclecost[i][j]);
        			handler.endArray();
        		}
        		handler.endArray();
        		handler.endElement();

        		handler.startElement("P");
        		handler.startArray();
        		for (int i = 0 ; i<P.length; i++) {
        			handler.startArray();
        			for (int j = 0 ; j<P[i].length; j++) 
        				handler.addIntItem(P[i][j]);
        			handler.endArray();
        		}
        		handler.endArray();
        		handler.endElement();

        		handler.startElement("expDemand");
        		handler.startArray();
        		for (int j = 0 ; j<expDemand.length ; j++)
        			handler.addNumItem(expDemand[j]);
        		handler.endArray();
        		handler.endElement();
        	}
        };

	}

