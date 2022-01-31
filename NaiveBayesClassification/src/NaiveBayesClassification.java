import java.io.*;
import java.util.*;

class training{

	static double training_matrix[] [] = new double[269][30];
	static int[] gene_index = new int[269];
	static String[] gene_class = new String[269];
	static int row=0;
	static int column=0;
	static double[] bmean = new double[30];
	static double[] mmean = new double[30];
	static double[] bvar = new double[30];
	static double[] mvar = new double[30];
	static double mprob=0;
	static double bprob=0;
	static int mcount =0;
	static int bcount =0;

	static void processtraining(String FileName)
	{
		try {

			//COLLECTING DATA FROM THE TRAINING SET
			BufferedReader readfile = new BufferedReader(new FileReader(FileName));
			String read;
			while ((read = readfile.readLine()) != null) 
			{
				String tokens[]=read.split("\\t");
				column = 0;
				for(int i=0;i<=1;i++)
				{
					String geneid=tokens[0];
					String geneclass=tokens[1];
					int gid=Integer.parseInt(geneid);
					gene_index[row] = gid;
					gene_class[row] = geneclass;
				}

				for(int i=2;i<tokens.length;i++)
				{
					double value=Double.parseDouble(tokens[i]);
					//System.out.println(value);
					training_matrix[row][column]=value;
					column++;

				}
				row++;
			}

			//CALCULATING CLASS PROBABILITY VALUES 

			/*** for class M ***/

			for(int i=0;i<269;i++)
			{
				if(gene_class[i].equalsIgnoreCase("M"))
				{
					mcount++;
				}

			}

			mprob = mcount/269.0;
			System.out.println(" Count of M = "+mcount);
			System.out.println(" Probability of M = "+mprob);

			/*** for class B ***/

			for(int i=0;i<269;i++)
			{
				if(gene_class[i].equalsIgnoreCase("B"))
				{
					bcount++;
				}

			}

			bprob = (bcount/269.0);
			System.out.println(" Count of B = "+bcount);
			System.out.println(" Probability of B = "+bprob);

			//CALCULATIING MEAN FOR EACH ATTRIBUTE IN EACH CLASS 

			ArrayList<Double> mlist = new ArrayList<Double>();
			ArrayList<Double> blist = new ArrayList<Double>();

			System.out.println(" Columns = "+column);
			System.out.println(" Rows = "+row);
			for(int i=0;i<column;i++)
			{
				blist.clear();
				mlist.clear();
				for (int j=0;j<row;j++)
				{
					if(gene_class[j].equalsIgnoreCase("M"))
					{

						mlist.add(training_matrix[j][i]);

					}
					if(gene_class[j].equalsIgnoreCase("B"))
					{
						blist.add(training_matrix[j][i]);

					}

				}	
				bmean[i] = mean(blist);
				mmean[i] = mean(mlist);
				bvar[i] = variance(blist);
				mvar[i] = variance(mlist);

			}

			readfile.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	public static double mean(ArrayList<Double> arr)
	{
		double sum=0,mean=0;
		for(int i=0;i<arr.size();i++)
		{

			sum +=arr.get(i);
		}
		mean = sum/arr.size();
		return mean;
	}

	public static double variance(ArrayList<Double> arr)
	{
		double sum=0,mean=0, var=0,sum_mean=0;
		for(int i=0;i<arr.size();i++)
		{
			sum +=arr.get(i);
		}
		mean = sum/arr.size();

		for(int i=0;i<arr.size();i++)
		{
			sum_mean= sum_mean + ((arr.get(i)-mean)*(arr.get(i)-mean));
		}

		var = sum_mean/arr.size();
		return var;
	}


}

class test{

	static int[] gene_index = new int[300];
	static String[] gene_class = new String[300];
	static double test_matrix[] [] = new double[300][30];

	static double mattprob[] = new double[30];
	static double battprob[] = new double[30];

	static double bprob = 0;
	static double mprob = 0;

	static int mcount =0;
	static int bcount =0;

	static int row =0;
	static int column=0;

	static double totmprob = 1;
	static double totbprob = 1;
	static double finalmprob = 0;
	static double finalbprob = 0;


	static void processtest(String FileName)
	{
		//COLLECTING DATA FROM THE TRAINING SET
		BufferedReader readfile;
		try {
			readfile = new BufferedReader(new FileReader(FileName));

			String read;
			while ((read = readfile.readLine()) != null) 
			{
				//read1=readfile2.readLine();
				String tokens[]=read.split("\\t");
				column = 0;
				for(int i=0;i<1;i++)
				{
					String geneid=tokens[0];
					int gid=Integer.parseInt(geneid);
					gene_index[row] = gid;
				}

				for(int i=1;i<tokens.length;i++)
				{
					double value=Double.parseDouble(tokens[i]);
					//System.out.println(value);
					test_matrix[row][column]=value;
					column++;

				}
				row++;
			}


			//CALCULATING THE PROBABILITY FOR EACH ATTRIBUTE


			double value = 0;

			for(int i=0;i<row;i++)
			{
				for(int j=0;j<column;j++)
				{

					value = test_matrix [i][j];
					mprob = mgauss(value,j);
					bprob = bgauss(value,j);

					mattprob[j] = mprob;
					battprob[j] = bprob;

				}

				finalmprob = mproduct(mattprob);
				finalbprob = bproduct(battprob);

				if(finalbprob > finalmprob)
				{
					gene_class[i] = "B";
					bcount++;
				}
				else if (finalbprob < finalmprob)
				{
					gene_class[i] = "M";
					mcount++;
				}

			}
			FileWriter outfile = new FileWriter("resultfile.txt");
			PrintWriter out = new PrintWriter(outfile);
			for(int i=0;i<300;i++)
			{
				System.out.println(i+1+".Gene no."+gene_index[i]+"====>"+gene_class[i]);
				out.println(gene_index[i]+"\t"+gene_class[i]);
			}
			out.close();

			System.out.println("The number of B genes = "+bcount);
			System.out.println("The number of M genes = "+mcount);

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}



	public static double bgauss(double v,int attr)
	{
		double mu = training.bmean[attr];
		double sigma2 = training.bvar[attr];
		double pi = Math.PI;

		double denom = Math.sqrt((2*pi*sigma2));
		double expo = (-(Math.pow((v-mu),2.0)))/(2*sigma2);

		double prob = (1/denom)*(Math.exp(expo));

		return prob;
	}

	public static double mgauss(double v,int attr)
	{
		double mu = training.mmean[attr];
		double sigma2 = training.mvar[attr];
		double pi = Math.PI;

		double denom = Math.sqrt((2*pi*sigma2));
		double expo = (-(Math.pow((v-mu),2.0)))/(2*sigma2);

		double prob = (1/denom)*(Math.exp(expo));

		return prob;
	}


	public static double mproduct(double[] arr)
	{
		double product=1;
		double totproduct =0;

		for(int i=0;i<arr.length;i++)
		{
			product = product*arr[i];
		}

		totproduct = product*training.mprob;
		return totproduct;
	}


	public static double bproduct(double[] arr)
	{
		double product=1;
		double totproduct =0;

		for(int i=0;i<arr.length;i++)
		{
			product = product*arr[i];
		}

		totproduct = product*training.bprob;
		return totproduct;
	}


}



public class NaiveBayesClassification {

	public static void main(String[] args)
	{
		System.out.println();
		System.out.println("TEST.TXT CLASSIFICATION");
		System.out.println("***********************");
		training.processtraining("training.txt");
		test.processtest("test.txt");

	}

}
