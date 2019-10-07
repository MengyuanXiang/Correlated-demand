package sample_opl;

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

import java.io.*;
import java.util.Arrays;

/**
 * To run from Mac OS
 * 
 * -Djava.library.path=/Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * Environment variable
 * 
 * DYLD_LIBRARY_PATH /Applications/CPLEX_Studio128/opl/bin/x86-64_osx/
 * 
 * @author gwren
 *
 */

public class Model {
   
   int Nbmonths;
   double[] expDemand; 
   double ordercost;
   double holdingcost; 
   double penaltycost;
   double unitcost;
   double initialStock;
   double[][] covariance;
   int Nbpartitions;
   
   String instanceIdentifier;
   
   public Model(
         int Nbmonths, 
         double[] expDemand, 
         double ordercost, 
         double holdingcost, 
         double penaltycost,
         double unitcost,
         double initialStock,
         double[][] covariance,
         int Nbpartitions,
         String instanceIdentifier){
      this.Nbmonths = Nbmonths;
      this.expDemand = expDemand;
      this.ordercost = ordercost;
      this.holdingcost = holdingcost;
      this.penaltycost = penaltycost;
      this.unitcost = unitcost;
      this.initialStock = initialStock;
      this.covariance = covariance;
      this.Nbpartitions = Nbpartitions;
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
  

   
	public static double[][] calucluteCovariance(double [] ExpDemand, double cv, double [] rho){
		double [] stdDemand =new double [ExpDemand.length];
		for (int i = 0; i < ExpDemand.length; i ++) {
			stdDemand[i] = cv*ExpDemand[i];
		}
		
		double [][] covariance = new double [ExpDemand.length][ExpDemand.length];
		/*
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
		*/
		
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
		return covariance;
		
	}
   
   public static void main(String[] args){
      //int Nbmonths = 10;
      //double[] expDemand = {200,50,100,300,150,200,100,50,200,150}; 
      //double[] stdDemand = {60,15,30,90,45,60,30,15,60,45}; 
      int Nbmonths = 15;
      double[] expDemand = {5,8,24,39,16,29,51,39,75,69,26,20,32,11,19}; 
      double cv = 0.2; 
      double ordercost = 100;
      double holdingcost = 1;
      double penaltycost = 10;
      double unitcost = 0;
      double initialStock = 0;
      double[] rho = {0.5};
      int Nbpartitions = 10;
      
      double[][] covariance = calucluteCovariance(expDemand, cv, rho);
      
      double lb = Double.NaN;
      double ub = Double.NaN;
      try{
         Model model = new Model(
               Nbmonths, 
               expDemand, 
               ordercost, 
               holdingcost,
               penaltycost,
               unitcost, 
               initialStock, 
               covariance,
               Nbpartitions,
               null);
         lb = model.solve("rs_milp_piecewise_penalty_lb");
         
      }catch(IloException e){
         e.printStackTrace();
      }
      
      try{
         Model model = new Model(
               Nbmonths, 
               expDemand, 
               ordercost, 
               holdingcost, 
               penaltycost,
               unitcost, 
               initialStock, 
               covariance,
               Nbpartitions,
               null);
         ub = model.solve("rs_milp_piecewise_penalty_ub");
      }catch(IloException e){
         e.printStackTrace();
      }
      
      //System.out.println(lb);
      //System.out.println(ub);
   }
   


public double solve(String model_name) throws IloException{
        //IloOplFactory.setDebugMode(true);
        IloOplFactory oplF = new IloOplFactory();
        IloOplErrorHandler errHandler = oplF.createOplErrorHandler(System.out);
        IloCplex cplex = oplF.createCplex();
        IloOplModelSource modelSource=oplF.createOplModelSourceFromStream(getMILPModelStream(new File("./opl_models/backorders/"+model_name+".mod")),model_name);
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
            boolean[] R = new boolean[Nbmonths];
            double[] S = new double[Nbmonths];
            for(int i = 0; i < Nbmonths; i++){
               //s[i] = cplex.getValue(opl.getElement("sValue").asNumVarMap().get(1+i));
               R[i] = Math.round(cplex.getValue(opl.getElement("purchase").asIntVarMap().get(1+i))) == 1 ? true : false;
               S[i] = cplex.getValue(opl.getElement("stock").asNumVarMap().get(1+i))+expDemand[i];
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
            return S[0];
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
         
            handler.startElement("Nbmonths");
            handler.addIntItem(Nbmonths);
            handler.endElement();
            
            handler.startElement("expDemand");
            handler.startArray();
            for (int j = 0 ; j<expDemand.length ; j++)
                handler.addNumItem(expDemand[j]);
            handler.endArray();
            handler.endElement();
            

            handler.startElement("ordercost");
            handler.addNumItem(ordercost);
            handler.endElement();
            
            handler.startElement("holdingcost");
            handler.addNumItem(holdingcost);
            handler.endElement();
            
            handler.startElement("penaltycost");
            handler.addNumItem(penaltycost);
            handler.endElement();
            
            handler.startElement("unitcost");
            handler.addNumItem(unitcost);
            handler.endElement();
            
            handler.startElement("initialStock");
            handler.addNumItem(initialStock);
            handler.endElement();
            
            handler.startElement("covariance");
            handler.startArray();
            for (int i = 0 ; i<covariance.length; i++) {
            	handler.startArray();
               	for (int j = 0 ; j<covariance.length; j++) 
                   handler.addNumItem(covariance[i][j]);
                handler.endArray();
            }
            handler.endArray();
            handler.endElement();
            
            handler.startElement("Nbpartitions");
            handler.addIntItem(Nbpartitions);
            handler.endElement();
            
            double[] means = linearization.LinearizationParameters.getMeans(Nbpartitions);
            handler.startElement("means");
            handler.startArray();
            for (int j = 0 ; j<means.length ; j++)
                handler.addNumItem(means[j]);
            handler.endArray();
            handler.endElement();
            
            double[] probabilities = linearization.LinearizationParameters.getProbabilities(Nbpartitions);
            handler.startElement("probabilities");
            handler.startArray();
            for (int j = 0 ; j<probabilities.length ; j++)
                handler.addNumItem(probabilities[j]);
            handler.endArray();
            handler.endElement();
            
            double error = linearization.LinearizationParameters.getError(Nbpartitions);
            handler.startElement("error");
            handler.addNumItem(error);
            handler.endElement();
        }
    };
}

