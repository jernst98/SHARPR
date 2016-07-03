package ernst.SHARPR;

import org.apache.commons.math3.stat.inference.*;
import org.apache.commons.math3.distribution.*;
import org.apache.commons.math3.linear.*;
import java.io.*;
import java.util.*;
import java.text.*;



public class SHARPR
{
    //SHARPR (Systematic High-resolution Activation and Repression Profiling using Reporters)
    //code written by Jason Ernst

    static int DEFAULT_THRESHOLD = 20;// 40;

    //static int DEFAULT_PSEUDO = 1;//10;


    static int DEFAULT_RNAPSEUDO = 1;
    static int DEFAULT_DNAPSEUDO =1;

    static double DEFAULT_VARPRIOR1 = 1;
    static double DEFAULT_VARPRIOR2 = 50;
    static double DEFAULT_VALUEBINS = 0.5;
    static int DEFAULT_MINCOUNTEXTREME = 500;

    static int DEFAULT_TILEEXTENSION = 10;
    static int DEFAULT_NUMOVERLAPBINS =5000;// 10000;
    static int DEFAULT_NUMOVERLAPBINS_MAX = 200;
    //static int DEFAULT_NORMALIZEWIDTH = 5000;

    String szrnafile;

    String szdnafile;

    String sznormoutfile;
    String sznormoutfilenoavg;

    int ndnacutoff;


    int ndnapseudocount;// = 1;

    int nrnapseudocount;

    String szconvertinputfile;

    String szconvertoutputfile;

    String szcombinefileset1;

    String szcombinefileset2;

    String szcombineoutputfile;

    int ntileindex;

    String szinterpolateinputfile;

    String szinterpolateoutputfile;
    String szdeconvolveinputfile;
    String szdeconvolveoutputfile;
    double dpriorvar;
    int ntilelength;
    int nstepsize;
    int numtiles; 

    boolean bstandardize;

    //int nnormalizewidth;

    String szdatafiles;

    boolean bmaxenrich;

    boolean brecentermean;
    boolean brecenterends;

    double dvarprior1;
    double dvarprior2;
    String szrnafilelist;
    String szdnafilelist;

    String szadjacentoutfile;


    /**
     * SHARPR - ExecuteAll
     */
    public SHARPR(String szrnafilelist, String szdnafilelist, int ntileindex, double dvarprior1, double dvarprior2, int ntilelength, 
                     int nstepsize, int numtiles, String szinterpolateoutputfile,
		  int ndnacutoff,int nrnapseudocount,int ndnapseudocount) throws IOException
    {
	this.szrnafilelist = szrnafilelist;
	this.szdnafilelist = szdnafilelist;
	this.ntileindex = ntileindex;
	this.dvarprior1 = dvarprior1;
	this.dvarprior2 = dvarprior2;
	this.ntilelength = ntilelength;
	this.nstepsize = nstepsize;
	this.numtiles = numtiles;
	this.szinterpolateoutputfile = szinterpolateoutputfile;
	this.ndnacutoff = ndnacutoff;
	this.nrnapseudocount = nrnapseudocount;
	this.ndnapseudocount = ndnapseudocount;
	this.bstandardize = true;

	executeAll();
    }

    /**
     * SHARPR - deconvolve
     */
    public SHARPR(String szdeconvolveinputfile, String szdeconvolveoutputfile, double dpriorvar, int ntilelength, 
                     int nstepsize, int numtiles, boolean bstandardize) throws IOException
    {
	this.szdeconvolveinputfile = szdeconvolveinputfile;
	this.szdeconvolveoutputfile = szdeconvolveoutputfile;
	this.dpriorvar = dpriorvar;
	this.ntilelength = ntilelength;
	this.nstepsize = nstepsize;
	this.numtiles = numtiles;

	this.bstandardize = bstandardize;

        deconvolveExact();
	//temporarily replacing with exact
	//deconvolve();
    }

    /**
     * SHARPR - normalize
     */
    public SHARPR(String szrnafile, String szdnafile, String sznormoutfile, String sznormoutfilenoavg, int ndnacutoff, int nrnapseudocount,int ndnapseudocount) throws IOException 
   {
       //this.nnormalizewidth = nnormalizewidth;
       this.szrnafile = szrnafile;
       this.szdnafile = szdnafile;
       this.sznormoutfile = sznormoutfile;
       this.sznormoutfilenoavg = sznormoutfilenoavg;
       this.ndnacutoff = ndnacutoff;
       this.nrnapseudocount = nrnapseudocount;
       this.ndnapseudocount = ndnapseudocount;

       normalize();
   }

    /**
     * SHARPR - converttable
     */
    public SHARPR(String szconvertinputfile, String szconvertoutputfile, int ntileindex,  boolean brecentermean, boolean brecenterends) throws IOException
    {
	this.szconvertinputfile = szconvertinputfile;
	this.szconvertoutputfile = szconvertoutputfile;
	this.ntileindex = ntileindex;
	this.brecentermean = brecentermean;
	this.brecenterends = brecenterends;

	converttable();
    }

    /**
     * SHARPR - linear interpolate
     */
    public SHARPR(int nstepsize, String szinterpolateinputfile, String szinterpolateoutputfile) throws IOException
    {
        this.szinterpolateinputfile = szinterpolateinputfile;
        this.szinterpolateoutputfile = szinterpolateoutputfile;
        this.nstepsize = nstepsize;

	linearInterpolate();
    }


    /**
     * SHARPR - combine files
     * szcombinefileset1 comma delimited files with one prior
     * szcombinefileset2 comma delimited file with prior 2
     */
    public SHARPR(String szcombinefileset1, String szcombinefileset2, String szcombineoutputfile) throws IOException 
    {
	this.szcombinefileset1 = szcombinefileset1;
	this.szcombinefileset2 = szcombinefileset2;
	this.szcombineoutputfile = szcombineoutputfile;

	combinefiles();
    }


    int nbasetooverlap;
    int numoverlapbins;
    String szreportercoordinates;
    String szoverlapcoordinates;
    String szbasepredictions;
    String szoverlapoutput;
    String szfeatureoutput;
    String szsummaryoutput;

    String szsequences;
    int nextension;
    double dvaluebins;
    boolean bvaluebins;
    int nmincountextreme;

    /**
     * SHARPR - enrichment
     *
     */
    public SHARPR(int nbasetooverlap, int numoverlapbins, String szreportercoordinates, String szoverlapcoordinates, String szbasepredictions, String szoverlapoutput,
                     String szfeatureoutput, boolean bmaxenrich, String szsummaryoutput, double dvaluebins, boolean bvaluebins, int nmincountextreme) throws IOException
    {
	this.nbasetooverlap = nbasetooverlap;
	this.numoverlapbins = numoverlapbins;
	this.szreportercoordinates = szreportercoordinates;
	this.szoverlapcoordinates = szoverlapcoordinates;
	this.szbasepredictions = szbasepredictions;
	this.szoverlapoutput = szoverlapoutput;
	this.szfeatureoutput = szfeatureoutput;
	this.bmaxenrich = bmaxenrich;
	this.szsummaryoutput = szsummaryoutput;
	this.dvaluebins = dvaluebins;
	this.bvaluebins = bvaluebins;
	this.nmincountextreme = nmincountextreme;

	computeenrichment();
       
	//computeenrichmentVALS();
    }



    /**
     *
     * SHARPR - adjacentsig
     */
    public SHARPR(String szdatafiles, int nstepsize, int numtiles, String szsequences, int nextension, 
                     String szreportercoordinates, String szadjacentoutfile) throws IOException
    {
	//                   new SHARPR(szdatafiles, numtiles,szsequences,nextension);
	this.szreportercoordinates = szreportercoordinates;
	this.szdatafiles = szdatafiles;
	this.nstepsize = nstepsize;
	this.numtiles = numtiles;
	this.szsequences = szsequences;
	this.nextension = nextension;
	this.szadjacentoutfile = szadjacentoutfile;
	adjacentchanges();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /*
   public static class RecRatio 
   { 
      double ddnacount;
      int nindex;
      double dratio;
      double drandom;

      RecRatio(double ddnacount, int nindex, double dratio, double drandom)
      {
	 this.ddnacount = ddnacount;
	 this.nindex = nindex;
	 this.dratio = dratio;
	 this.drandom = drandom;
      }
   }


   public static class RecRatioCompare implements Comparator
   {
      public int compare(Object o1, Object o2)
      {
         RecRatio r1 = (RecRatio) o1;
	 RecRatio r2 = (RecRatio) o2;

	 if (r1.ddnacount < r2.ddnacount)
	 {
	    return -1;
	 }
	 else if (r1.ddnacount > r2.ddnacount)
	 {
	    return 1;
	 }
	 else if (r1.drandom < r2.drandom)
	 {
	    return -1;
	 }
	 else if (r1.drandom > r2.drandom) 
         {
	    return 1;
	 }
	 else
	 {
	    return 0;
	 }
      }
   }
    */

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void executeAll() throws IOException
    {
       if (ntilelength % nstepsize != 0)
	   throw new IllegalArgumentException("tilelength is not divisible by stepsize values are "+ntilelength+" "+nstepsize);


       double dlog2 = Math.log(2);


       ////////////////////////////////////////////////////////////////////////////////////////
       ///////////////////////////////////////////////////////////////////////////////////////
       ////////// normalize part

       StringTokenizer strnafiles = new StringTokenizer(szrnafilelist,",");
       StringTokenizer stdnafiles = new StringTokenizer(szdnafilelist,",");
       if (strnafiles.countTokens() != stdnafiles.countTokens())
       {
	   throw new IllegalArgumentException("Mismatch in the number of RNA and DNA files specified "+strnafiles.countTokens()+" vs. "+stdnafiles.countTokens());
       }
       int numfiles = strnafiles.countTokens();
       String[][] normratio = null;

       String[][] szreporterIDA = null;

       for (int nfile = 0; nfile < numfiles; nfile++)
       {
	   String szdnafile = (String) stdnafiles.nextToken();
	   String szrnafile = (String) strnafiles.nextToken();

          BufferedReader brdna = Util.getBufferedReader(szdnafile);

          String szLine_rna;
          String szLine_dna = brdna.readLine();

          //String szLine;

          int numlines = 0;
	  if (szLine_dna == null)
	  {
	      throw new IllegalArgumentException(szdnafile+" is empty!");
	  }
          StringTokenizer st = new StringTokenizer(szLine_dna,"\t");
           //the number of columns corresponds to the number of headers - assumes no header for id column
          int numcols = st.countTokens();
          double dsumdna = 0;
          double dsumrna = 0;


           //storing the probes and the total dna count
          //ArrayList alprobes = new ArrayList();
          while ((szLine_dna = brdna.readLine())!=null)
          {
              st = new StringTokenizer(szLine_dna,"\t");
              String sztoken = st.nextToken();
              //alprobes.add(sztoken);
              while (st.hasMoreTokens())
              {
                 dsumdna += Double.parseDouble(st.nextToken())+ndnapseudocount;
              }
              numlines++;
          }
          brdna.close();

          brdna = Util.getBufferedReader(szdnafile);
          brdna.readLine();  //flushing the header

          BufferedReader brrna = Util.getBufferedReader(szrnafile);

          brrna.readLine();  //flushing the header
          int nline = 0;

          int[][] dnacount = new int[numlines][numcols]; //stores the original count value
          int[][] rnacount = new int[numlines][numcols];

          double[][] dnaval = new double[numlines][numcols]; //stores the original count value
          double[][] rnaval = new double[numlines][numcols];
  
          double[][] ratio = new double[numlines][numcols];

	  if (nfile == 0)
	  {
	     szreporterIDA = new String[numfiles][numlines];
             normratio = new String[numfiles][numlines];
	  }

          while ((szLine_rna = brrna.readLine())!=null)
          {
             szLine_dna = brdna.readLine();

	     StringTokenizer strna = new StringTokenizer(szLine_rna,"\t");
	     String szrnaid = strna.nextToken();

	     if (szLine_dna == null)
	     {
		 throw new IllegalArgumentException("There are more lines in "+szrnafile+" than in  "+szdnafile+" expecting the same number");
	     }
 	     StringTokenizer stdna = new StringTokenizer(szLine_dna,"\t");
	     String szdnaid = stdna.nextToken();

	     if (!szrnaid.equals(szdnaid))
	     {
	        throw new IllegalArgumentException("RNA file "+szrnafile+" does not match DNA file "+szdnafile+" with ID "+szrnaid+" "+szdnaid);
	     }
	     szreporterIDA[nfile][nline] = szrnaid;

	     int ncol = 0;
	     while (strna.hasMoreTokens())
	     {
	        int nrnaval = nrnapseudocount+Integer.parseInt(strna.nextToken());
	        int ndnaval = ndnapseudocount+Integer.parseInt(stdna.nextToken());
	        dsumrna += nrnaval;
	        rnacount[nline][ncol] = nrnaval;
	        dnacount[nline][ncol] = ndnaval;
	        ncol++;
	     }
	     nline++;
	  }
	  brrna.close();
	  brdna.close();
      

          for (nline = 0; nline < numlines; nline++)
          {
             for (int ncol = 0; ncol < numcols; ncol++)
	     {
	        rnaval[nline][ncol] = Math.log(rnacount[nline][ncol]/(double) dsumrna)/dlog2;
		dnaval[nline][ncol] = Math.log(dnacount[nline][ncol]/(double) dsumdna)/dlog2;
		ratio[nline][ncol] = rnaval[nline][ncol] - dnaval[nline][ncol];
	     }

	     ArrayList ratiosAL = new ArrayList();
	     for (int ncol = 0; ncol < numcols; ncol++)
	     {
	        if ((dnacount[nline][ncol]>=(ndnacutoff+ndnapseudocount)))
	        {
	           double dratioval = ratio[nline][ncol];
	           ratiosAL.add(Double.valueOf(dratioval));
		}
	     }

	     //copies ratio values into array so can easily be sorted
             double[] copyratios = new double[ratiosAL.size()];

	     for (int nj = 0; nj < copyratios.length; nj++)
	     {
                copyratios[nj] = ((Double) ratiosAL.get(nj)).doubleValue();
	     }
	     Arrays.sort(copyratios);

	     if (copyratios.length == 0)
	     {
		normratio[nfile][nline] = null;
	     }
	     else if (copyratios.length % 2 == 0)
	     {
		normratio[nfile][nline] =""+ ((copyratios[copyratios.length/2-1]+copyratios[copyratios.length/2])/2.0);
                //even number of ratio, takes the averge since in log space
	     }
	     else
	     {
                normratio[nfile][nline] = ""+copyratios[copyratios.length/2];
	     }
	  }
       }

       ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       ///now doing the conversion into table format

       String[][][] convertedvals = new String[numfiles][][];
       String[][] convertedElementIDs = new String[numfiles][];

       for (int nfile = 0; nfile < numfiles; nfile++)
       {
          String szLine;

          int ni = 0;
          boolean bfirst = true;

          int nprevpos=-1;
          String szprevlabel="";

          HashMap hmlabels = new HashMap();
          ArrayList allabels = new ArrayList();
          int nlabelindex = 0;
          int nmaxpos = -1;

	  for (int nreporter = 0; nreporter < szreporterIDA[nfile].length; nreporter++)
          {
	     StringTokenizer stu = new StringTokenizer(szreporterIDA[nfile][nreporter],"_");

	     int npos= -1;
	     int nuindex = 0;
  	     boolean bufirst = true;
	     String szlabel = "";
	     while (stu.hasMoreTokens())
	     {
	        if (nuindex == ntileindex)
	        {
	           npos = Integer.parseInt(stu.nextToken());
		   if (npos > nmaxpos)
		   {
		      nmaxpos = npos;
		   }
	        }
	        else
	        {
	           if (bufirst)
	           {
		      szlabel = stu.nextToken();
		      bufirst = false;
	           }
	           else
	           {
		     // if (nuindex <= ntileindex)
		      szlabel += "_" + stu.nextToken();
		    // else 
		    //  stu.nextToken();
		   }
	         }
	         nuindex++;
	     }

	     Integer objcount = (Integer) hmlabels.get(szlabel);
	     if (objcount == null)
	     {
	        allabels.add(szlabel);
	        hmlabels.put(szlabel, Integer.valueOf(nlabelindex));
	        nlabelindex++;
	     }
	  }
          //br.close();


          convertedvals[nfile] = new String[nlabelindex][nmaxpos+1];
	  convertedElementIDs[nfile] = new String[nlabelindex];
	  for (int na = 0; na < allabels.size(); na++)
	  {
	      convertedElementIDs[nfile][na] = (String) allabels.get(na);
	  }

	  for (int nreporter = 0; nreporter < normratio[nfile].length; nreporter++)
          {
	     String szID = szreporterIDA[nfile][nreporter];
	     StringTokenizer stu = new StringTokenizer(szID,"_");
 
 	     int npos= -1;
	     int nuindex = 0;
	     boolean bufirst = true;
	     String szlabel = "";
	     while (stu.hasMoreTokens())
	     {
	        if (nuindex == ntileindex)
	        {
	           npos = Integer.parseInt(stu.nextToken());
	        }
	        else
	        {
	           if (bufirst)
	           {
		      szlabel = stu.nextToken();
		      bufirst = false;
	           }
	           else
	           {
		     // if (nuindex <= ntileindex)
		      szlabel += "_" + stu.nextToken();
		    // else 
		    //  stu.nextToken();
		   }
		}
	        nuindex++;
	     }


	     int nrowindex = ((Integer) hmlabels.get(szlabel)).intValue();
             //String szval = st.nextToken();
             //int ncount = Integer.parseInt(st.nextToken());
	     if (normratio[nfile][nreporter] != null)
	     {
		 convertedvals[nfile][nrowindex][npos] = normratio[nfile][nreporter];// szval;
	     }
	  }
       }


       /////////////////////////////////////////////////////////////////////////////////
       /////////////////////////////////////////////////////////////////////////////////
       ////////////// infer part

       NumberFormat nf  = NumberFormat.getInstance();
       nf.setMaximumFractionDigits(3);

 
       int numhiddenvalstotal = ((ntilelength) + nstepsize * (numtiles-1))/nstepsize;
       int numhiddenvalstile = ntilelength/nstepsize;

       double[] priorvalsA = new double[2];
       priorvalsA[0] = dvarprior1;
       priorvalsA[1] = dvarprior2;
       String[][][][] infervalsA = new String[priorvalsA.length][numfiles][][];

       for (int nfile = 0; nfile < numfiles; nfile++)
       {
	   String szLine;

	   int numlines = convertedvals[nfile].length;

	   double[][] data = new double[numlines][numtiles];
	   boolean[][] bpresent = new boolean[numlines][numtiles];

	   for (int npriorindex = 0; npriorindex < infervalsA.length; npriorindex++)
	   {
              infervalsA[npriorindex][nfile] = new String[numlines][numhiddenvalstotal];
	   }


	   //ArrayList alvals = new ArrayList();
	   int npos = 0;
	   double dval = 0;
	   double dvalsq = 0;
 
	   for (int nrow = 0; nrow < numlines; nrow++)
           {
	       for (int ncol = 0; ncol <  convertedvals[nfile][nrow].length; ncol++)
	       {
	          if (ncol >= numtiles)
	          {
		      throw new IllegalArgumentException("We have more entries for "+convertedElementIDs[nfile][nrow]+" "+(ncol+1)+" than number of tile positions specified "+numtiles);
		  }

		  if (convertedvals[nfile][nrow][ncol] != null)
	          {
		     data[nrow][ncol] = Double.parseDouble(convertedvals[nfile][nrow][ncol]);//sztoken);
		     bpresent[nrow][ncol] = true;

		     npos++;
		     dval += data[nrow][ncol];
		     dvalsq += data[nrow][ncol]*data[nrow][ncol];
		  }
	       }
	   }
   
	   double dinitmean = dval/npos;
	   double dnormvar = dvalsq/npos -(dval/npos)*(dval/npos);

	   /*
           boolean BZSCORE = false;//true;

           if (BZSCORE)
           {
	      double dstd = Math.sqrt(dnormvar);

	      for (int nrow = 0; nrow < data.length; nrow++)
	      {
	         for (int ncol = 0; ncol < data[nrow].length; ncol++)
	         {
		    if (bpresent[nrow][ncol])
		    {
		       data[nrow][ncol] = (data[nrow][ncol]-dinitmean)/dstd;
		    }
		 }
	      }
	      dinitmean = 0;
	      dnormvar = 1;
	   }
	   */

	   for (int npriorindex = 0; npriorindex < priorvalsA.length; npriorindex++)
	   {
	       //double[] datavals = new double[numhiddenvalstotal];

	       //double[] datavalssum = new double[numhiddenvalstotal];
	       double[][] dhiddenvalA = new double[numlines][numhiddenvalstotal];

	       double[][] meanX = new double[numhiddenvalstotal][1];

               for (int nj = 0; nj < numhiddenvalstotal; nj++)
               {
	          meanX[nj][0] = dinitmean;
	       }

	       RealMatrix theRealMatrixMeanX = MatrixUtils.createRealMatrix(meanX);

               for (int nrow = 0; nrow < data.length; nrow++)
               {
		   //for (int nval = 0; nval < datavals.length; nval++)
		   //{
	           //  datavals[nval] = dinitmean;
		   //   datavalssum[nval] = 0;
		   //}
	          double[] data_nrow = data[nrow];
	          boolean[] bpresent_nrow = bpresent[nrow];

	          int numpresenttiles = 0;
	          for (int nk = 0; nk < bpresent_nrow.length; nk++)
	          {
	             if (bpresent_nrow[nk])
	             {
		        numpresenttiles++;
		     }
		  }

	         if (numpresenttiles == 0)
	         {
	            for (int nk = 0; nk < dhiddenvalA[nrow].length; nk++)
	            {
	               dhiddenvalA[nrow][nk] = dinitmean;
		    }
		 }
	         else
	         {
		    double[][] covarXY = new double[numhiddenvalstotal][numpresenttiles];
		    double[][] covarYY = new double[numpresenttiles][numpresenttiles];
		    double[][] observedY = new double[numpresenttiles][1];
		    double[][] meanY = new double[numpresenttiles][1];

		    int[] presentindiciesA = new int[numpresenttiles];

		    int npresenttile = 0;
	            for (int ntile = 0; ntile < numtiles; ntile++)
	            {
	               if (bpresent_nrow[ntile])
	               {
		          observedY[npresenttile][0] = data_nrow[ntile];
		  	  meanY[npresenttile][0] = dinitmean;

			  presentindiciesA[npresenttile] = ntile;

	  	          for (int nxval = ntile; nxval < ntile + numhiddenvalstile; nxval++)
		          {
	                     covarXY[nxval][npresenttile] = priorvalsA[npriorindex]/(double) numhiddenvalstile;		
			  }
			  npresenttile++;
		       }
		    }

 	            for (npresenttile = 0; npresenttile < numpresenttiles; npresenttile++)
	            {
	               for (int npresenttile2 = 0; npresenttile2 < numpresenttiles; npresenttile2++)
	               {
		          if (npresenttile == npresenttile2)
		          {
		             covarYY[npresenttile][npresenttile2] = dnormvar + priorvalsA[npriorindex]/(double) numhiddenvalstile;
			  }
		          else
		          {
		             if (Math.abs(presentindiciesA[npresenttile]-presentindiciesA [npresenttile2]) < numhiddenvalstile)
		             {
			        covarYY[npresenttile][npresenttile2] = priorvalsA[npriorindex] * 
                                                             (numhiddenvalstile-Math.abs(presentindiciesA[npresenttile]-presentindiciesA[npresenttile2]))/
                                                                                     (double) (numhiddenvalstile*numhiddenvalstile);

			     }
			  }
		       }
		    } 
	      
		    RealMatrix theRealMatrixMeanY = MatrixUtils.createRealMatrix(meanY);
		    RealMatrix theRealMatrixObservedY = MatrixUtils.createRealMatrix(observedY);
		    RealMatrix theRealMatrixCovarXY = MatrixUtils.createRealMatrix(covarXY);
		    RealMatrix theRealMatrixCovarYY = MatrixUtils.createRealMatrix(covarYY);

                    dhiddenvalA[nrow] = theRealMatrixMeanX.add(theRealMatrixCovarXY.multiply(MatrixUtils.inverse(theRealMatrixCovarYY).multiply(							                                                    theRealMatrixObservedY.subtract(theRealMatrixMeanY)))).getColumn(0);

		 }
	       }

               if (bstandardize)
               {
                  double dsumhidden = 0;
                  double dsumhiddensq = 0;
	          int numels = 0;
	          for (int nrowindex = 0; nrowindex < dhiddenvalA.length; nrowindex++)
	          {
	             for (int nk = 0; nk < dhiddenvalA[nrowindex].length; nk++)
	             {
		        double dcurrval = dhiddenvalA[nrowindex][nk];

		        dsumhidden += dcurrval;
		        dsumhiddensq += dcurrval * dcurrval;
		        numels++;
		     }
		  }

   	          double dmeanhidden = dsumhidden/numels;
	          double dvarhidden = dsumhiddensq/numels - dmeanhidden * dmeanhidden;
	          double dstdhidden = Math.sqrt(dvarhidden);

                  for (int nrowindex = 0; nrowindex < dhiddenvalA.length; nrowindex++)
                  {
		     double[] dhiddenvalA_nrowindex = dhiddenvalA[nrowindex];
	             for (int nk = 0; nk < dhiddenvalA_nrowindex.length; nk++)
                     {
		        infervalsA[npriorindex][nfile][nrowindex][nk] = nf.format((dhiddenvalA_nrowindex[nk]-dmeanhidden)/dstdhidden);
       		     }
		  }
	       }
               else
               {
                  for (int nrowindex = 0; nrowindex < dhiddenvalA.length; nrowindex++)
                  {
		      double[] dhiddenvalA_nrowindex = dhiddenvalA[nrowindex];
	             for (int nk = 0; nk < dhiddenvalA_nrowindex.length; nk++)
                     {
		        infervalsA[npriorindex][nfile][nrowindex][nk] = nf.format(dhiddenvalA_nrowindex[nk]);
		     }
		  }
	       }
	   }
       }


       /////////////////////////////////////////////////////////////////////////////////
       /////////////////////////////////////////////////////////////////////////////////
       ////////////// combine part
       String[][] combinevalsA = new String[infervalsA[0][0].length][infervalsA[0][0][0].length];


	String szLine;
	for (int nrow = 0; nrow < combinevalsA.length; nrow++)
	{
	   String szID = convertedElementIDs[0][nrow];//(String) stA1[0].nextToken();
	   for (int nk = 1; nk < convertedElementIDs.length; nk++)
	   {
	      String szcurrID =  convertedElementIDs[nk][nrow];// stA1[nk].nextToken();
	      if (!szcurrID.equals(szID))
	      {
	         throw new IllegalArgumentException("We had different IDs "+szcurrID+" and "+szID+" on corresponding positions!");
	      }
	   }

	   for (int ncol = 0; ncol < combinevalsA[0].length; ncol++)
	   {
	       double dsum1 = Double.parseDouble(infervalsA[0][0][nrow][ncol]);//stA1[0].nextToken());

		for (int nk =1; nk < numfiles; nk++)
		{
                    dsum1 += Double.parseDouble(infervalsA[0][nk][nrow][ncol]);//stA1[nk].nextToken());
		}

		double davg1 = dsum1/numfiles;//stA1.length;

		double doutput = davg1;


		 double dsum2 = Double.parseDouble(infervalsA[1][0][nrow][ncol]);

	         for (int nk =1; nk < numfiles; nk++)
	         {
		    dsum2 += Double.parseDouble(infervalsA[1][nk][nrow][ncol]);// stA2[nk].nextToken());
	         }

		 double davg2 = dsum2/numfiles;//stA2.length;

	         if (doutput >= 0)
	         {
	            doutput = Math.max(Math.min(davg2, doutput), 0);
	         }
	         else 
	         {
	            doutput = Math.min(Math.max(davg2, doutput), 0);
		 }

		//}
		 combinevalsA[nrow][ncol] = nf.format(doutput);
	    }
	}

       ///////////////////////////////////////////////////////////////////////////////
       //////////////////////////////////////////////////////////////////////////////
       /////////////  interpolate part

	PrintWriter pw = new PrintWriter(new FileWriter(szinterpolateoutputfile));
	//String szLine;

	double[] data = new double[combinevalsA[0].length];// null;

	for (int nrow = 0; nrow < combinevalsA.length; nrow++)
	{

	   pw.print(convertedElementIDs[0][nrow]);// st.nextToken());

	   for (int nindex  = 0; nindex < combinevalsA[nrow].length; nindex++)
	   {
	       //System.out.println(nrow+"\t"+nindex+"\t"+combinevalsA[nrow][nindex]);
	       data[nindex] = Double.parseDouble(combinevalsA[nrow][nindex]);// st.nextToken());
	   }

	   int numpos = nstepsize*data.length;
  	   for (int npos = 0; npos < numpos; npos++)
	   {
	      double dsmoothval;

	      int nhiddenindex = npos/nstepsize;
	      int nhiddenpos = npos % nstepsize;

	      //4

              if (((nhiddenindex ==0)&&(nhiddenpos <nstepsize/2))||(nhiddenpos == nstepsize/2)||((nhiddenindex == data.length-1)&&(nhiddenpos>=nstepsize/2)))
	      {
	         dsmoothval = data[nhiddenindex];
	      }
	      else if (nhiddenpos < nstepsize/2)
	      {
		  //0,1
		  //0--> 3/4, 1/4
		  //1-->1, 0
	         dsmoothval = ((nstepsize/2+nhiddenpos+1)/(double) nstepsize)*data[nhiddenindex]+(1-(nstepsize/2+nhiddenpos+1)/(double) nstepsize)*data[nhiddenindex-1];
	      }
	      else
	      {
		  //2-->1,0
		  //3-->3/4,1/4
	         dsmoothval = (1-(nhiddenpos-nstepsize/2)/(double) nstepsize)*data[nhiddenindex]+((nhiddenpos-nstepsize/2)/(double) nstepsize)*data[nhiddenindex+1];
	      }

	      pw.print("\t"+nf.format(dsmoothval));
	   }
	   pw.println();
	}
	pw.close();

    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     *
     */
    public void computeenrichment() throws IOException
    {

	BufferedReader brreportercoordinates = Util.getBufferedReader(szreportercoordinates);
	BufferedReader broverlapcoordinates = Util.getBufferedReader(szoverlapcoordinates);

	PrintWriter pwfeature = null;
	if (szfeatureoutput != null)
	{
	    pwfeature = new PrintWriter(new FileWriter(szfeatureoutput));
	}

	HashMap hmreportercoord = new HashMap();
	HashMap hmfeature = new HashMap();
	String szLine;

	NumberFormat nf = NumberFormat.getInstance();
	nf.setMaximumFractionDigits(5);

	//storing the feature information

	while ((szLine = brreportercoordinates.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
	    String szid = st.nextToken();
	    String szchrom =  st.nextToken();
	    int nbegin = Integer.parseInt(st.nextToken());
	    int nend = Integer.parseInt(st.nextToken());

	    int npos = 0;
	    for (int ncoord = nbegin; ncoord < nend; ncoord++)
	    {
		if ((nbasetooverlap == -1)||(nbasetooverlap == npos))
		{
		    hmfeature.put(szchrom+"\t"+ncoord, Boolean.valueOf(false));
		}
		npos++;
	    }

	    hmreportercoord.put(szid,szchrom+"\t"+nbegin+"\t"+nend);
	}
	brreportercoordinates.close();

	while ((szLine = broverlapcoordinates.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
	    String szchrom = st.nextToken();
            int nbegin = Integer.parseInt(st.nextToken());
            int nend = Integer.parseInt(st.nextToken());

            for (int ncoord = nbegin; ncoord < nend; ncoord++)
	    {
	       if (hmfeature.get(szchrom+"\t"+ncoord)!=null)
	       {
                  hmfeature.put(szchrom+"\t"+ncoord, Boolean.valueOf(true));
	       }
	    }	     
	}
	broverlapcoordinates.close();

	BufferedReader brbasepredictions = Util.getBufferedReader(szbasepredictions);

	Random theRandom = new Random(3433);
	Random theRandom2 = new Random(456);
	ArrayList albase = new ArrayList();
	ArrayList albaseabs = new ArrayList();

	int ncenterhit = 0;
	int ncentertotal = 0;
        double dsumabspredictnohit = 0;
	double dsumabspredicthit = 0;
	double dsumabspredict = 0;
	double dsumcenterabspredicthit = 0;
        double dsumcenterabspredict = 0;


	//double dmaxabsbestpredict = 0;



	//boolean bmax = false; //variable we need to toggle, true if want max, false all positions


	while ((szLine = brbasepredictions.readLine())!=null)
	{

	    StringTokenizer st = new StringTokenizer(szLine,"\t");
	    String szid = st.nextToken();

	    String szcoords = (String) hmreportercoord.get(szid);
	    StringTokenizer stcoords = new StringTokenizer(szcoords,"\t");
	    String szchrom = stcoords.nextToken();
	    int nbegin = Integer.parseInt(stcoords.nextToken());
	    int nend = Integer.parseInt(stcoords.nextToken());

	    int npos = 0;

	    int ncenterpos = (nbegin+nend)/2-nbegin;

	    if (szfeatureoutput != null)
	    {
		pwfeature.print(szid);
	    }

	    double dmaxabsbestpredict = -1;
	    double dmaxabsbestpredictnonabsval = 0;
	    boolean bmaxhit = false;

	    for (int ncoord = nbegin; ncoord < nend; ncoord++)
	    {
		if ((nbasetooverlap == -1)||(nbasetooverlap == npos))
		{

		   boolean bhit = ((Boolean) hmfeature.get(szchrom+"\t"+ncoord)).booleanValue();

		   if (szfeatureoutput != null)
		   {
		       if (bhit)
		       {
			   pwfeature.print("\t1");
		       }
		       else
		       {
			   pwfeature.print("\t0");
		       }
		   }

		   double dpredictval = Double.parseDouble(st.nextToken());


		   //else 
		   //{
                   //    dsumabspredictnohit += Math.abs(dpredictval);
		   //}                       

		   if (bmaxenrich)
		   {
		       if (Math.abs(dpredictval) > dmaxabsbestpredict)
		       {
			   dmaxabsbestpredict = Math.abs(dpredictval);
			   dmaxabsbestpredictnonabsval = dpredictval;
			   bmaxhit = bhit;
		       }
		   }
		   else //if (!bmaxenrich)
		   {
		       if (bhit)
		       {
		          dsumabspredicthit += Math.abs(dpredictval);
		       }
		       dsumabspredict += Math.abs(dpredictval);

		      albase.add(new OverlapRec(dpredictval, bhit, theRandom.nextDouble()));
		      albaseabs.add(new OverlapRec(Math.abs(dpredictval), bhit, theRandom2.nextDouble()));
		   }

		   if (npos == ncenterpos)
	           {
		       if (bhit)
		       {
		          ncenterhit++;
		       
		          dsumcenterabspredicthit += Math.abs(dpredictval);
		       }
		       dsumcenterabspredict += Math.abs(dpredictval);
		       ncentertotal++;		     
		       //System.out.println("==>"+Math.abs(dpredictval)+"\t"+ncentertotal);
  
 		   }
		}
		else
		{
		    st.nextToken();
		}


		npos++;
	    }

	    if (bmaxenrich)
	    {
		albase.add(new OverlapRec(dmaxabsbestpredictnonabsval, bmaxhit, theRandom.nextDouble()));
		albaseabs.add(new OverlapRec(dmaxabsbestpredict, bmaxhit, theRandom2.nextDouble()));

		if (bmaxhit)
		{
		   dsumabspredicthit += dmaxabsbestpredict;
		}
		dsumabspredict += dmaxabsbestpredict;
	    }

	    if (szfeatureoutput != null)
	    {
		pwfeature.println();
	    }

	}
	brbasepredictions.close();

	if (szfeatureoutput !=null)
	{
	    pwfeature.close();
	}

	OverlapRec[] baseA = new OverlapRec[albase.size()];
	for (int nk = 0; nk < baseA.length; nk++)
	{
	    baseA[nk] = (OverlapRec) albase.get(nk);
	}
	Arrays.sort(baseA, new OverlapRecCompare());

        OverlapRec[] baseabsA = new OverlapRec[albaseabs.size()];
        for (int nk = 0; nk < baseabsA.length; nk++)
	{
      	   baseabsA[nk] = (OverlapRec) albaseabs.get(nk);
	}
        Arrays.sort(baseabsA, new OverlapRecCompare());


        double dbinsum = 0;
        int nbincount = 0;

	int nsumabove = 0;
	int nsumnotabove = 0;
	int ncurrabove = 0;
	int ncurrnotabove = 0;
	int nreversesumabove = 0;
	int nreversesumnotabove = 0;

        int ntotabove = 0;
        int ntotnotabove = 0;
        for (int nindex =0; nindex < baseA.length; nindex++)
	{
	    if (baseA[nindex].bhit)
	    {
		ntotabove++;
	    }
	    else
	    {
		ntotnotabove++;
	    }
	}


	int ntp = 0;
	int nfp = 0;
	int nlastfp = 0;
	double dauc = 0;
	double dauc005 = 0;
	double dauc01 = 0;
	double dauc05 = 0;

	//for (int nindex =0; nindex < baseabsA.length; nindex++)
	for (int nindex =0; nindex < baseA.length; nindex++)
	{
	    //if (baseabsA[nindex].bhit)
	    if (baseA[nindex].bhit)
	    {
		ntp++;
	    }
	    else
	    {
		nfp++;
	    }

	    double dincrement = ((nfp-nlastfp)/(double) (ntotnotabove))*(ntp/(double) ntotabove);
	    dauc += dincrement;

	    if (nfp/(double) ntotnotabove <= 0.05)
	    {
		dauc05 += dincrement;

		if (nfp/(double) ntotnotabove <= 0.01)
		{
		    dauc01 += dincrement;

		    if (nfp/(double) ntotnotabove <= 0.005)
		    {
			dauc005 += dincrement;
		    }

		}
	    }

	    nlastfp = nfp;
	}


        ntp = 0;
        nfp = 0;
        nlastfp = 0;
        double daucabs = 0;
        double dauc01abs = 0;
        double dauc05abs = 0;
	double dauc005abs = 0;

	for (int nindex =0; nindex < baseabsA.length; nindex++)
	{
	    if (baseabsA[nindex].bhit)
	    {
		ntp++;
	    }
	    else
	    {
		nfp++;
	    }

	    double dincrement = ((nfp-nlastfp)/(double) (ntotnotabove))*(ntp/(double) ntotabove);
	    daucabs += dincrement;

	    if (nfp/(double) ntotnotabove <= 0.05)
	    {
		dauc05abs += dincrement;

		if (nfp/(double) ntotnotabove <= 0.01)
		{
		    dauc01abs += dincrement;

                    if (nfp/(double) ntotnotabove <= 0.005)
		    {
		       dauc005abs += dincrement;
		    }
		}
	    }

	    nlastfp = nfp;
	}



        ntp = 0;
        nfp = 0;
        nlastfp = 0;
        double daucrev = 0;
        double dauc01rev = 0;
        double dauc05rev = 0;
	double dauc005rev = 0;

	for (int nindex = baseA.length-1; nindex >= 0; nindex--)
	{
	    //if (baseabsA[nindex].bhit)
	    if (baseA[nindex].bhit)
	    {
		ntp++;
	    }
	    else
	    {
		nfp++;
	    }

	    double dincrement = ((nfp-nlastfp)/(double) (ntotnotabove))*(ntp/(double) ntotabove);
	    daucrev += dincrement;

	    if (nfp/(double) ntotnotabove <= 0.05)
	    {
		dauc05rev += dincrement;

		if (nfp/(double) ntotnotabove <= 0.01)
		{
		    dauc01rev += dincrement;

                    if (nfp/(double) ntotnotabove <= 0.005)
		    {
		       dauc005rev += dincrement;
		    }
		}
	    }

	    nlastfp = nfp;
	}

	//System.out.println("STATS\t"+szbasepredictions+"\t"+szoverlapcoordinates+"\t"+(ntotabove/(double) (ntotabove+ntotnotabove))+"\t"+(dsumabspredicthit/ntotabove)+"\t"+
	//   (dsumabspredict/(double) (ntotabove+ntotnotabove))+"\t"+(dsumcenterabspredicthit/ncentertotal));

	PrintWriter pwstats = new PrintWriter(new FileWriter(szoverlapoutput));
        double dtotratio= ntotabove/(double) (ntotabove+ntotnotabove);
	double dcenterratio = ncenterhit/(double) ncentertotal;


	if (szsummaryoutput != null)
	{
	    PrintWriter pwsummary = new PrintWriter(new FileWriter(szsummaryoutput));
	    pwsummary.println("Base prediction file\t"+szbasepredictions);
	    pwsummary.println("Overlap coordinates file\t"+szoverlapcoordinates);
	    pwsummary.println("Activating Ranking AUC\t"+dauc);
	    pwsummary.println("Activating Ranking AUC 5%\t"+dauc05);
	    pwsummary.println("Activating Ranking AUC 1%\t"+dauc01);
            pwsummary.println("Activating Ranking AUC 0.5%\t"+dauc005);

            pwsummary.println("Repressive AUC\t"+daucrev);
            pwsummary.println("Repressive AUC 5%\t"+dauc05rev);
            pwsummary.println("Repressive AUC 1%\t"+dauc01rev);
            pwsummary.println("Repressive AUC 0.5%\t"+dauc005rev);

            pwsummary.println("Absolute AUC\t"+daucabs);
            pwsummary.println("Absolute AUC 5%\t"+dauc05abs);
            pwsummary.println("Absolute AUC 1%\t"+dauc01abs);
            pwsummary.println("Absolute AUC 0.5%\t"+dauc005abs);

            pwsummary.println("Fraction overlap overall\t"+dtotratio);
	    pwsummary.println("Average abs. activity of overlapped bases\t"+(dsumabspredicthit/(double) ntotabove));
            pwsummary.println("Average abs. activity overall\t"+(dsumabspredict/(double) (ntotabove+ntotnotabove)));
            pwsummary.println("Fraction overlap center\t"+dcenterratio);
            pwsummary.println("Average abs. activity of overlapped bases center\t"+(dsumcenterabspredicthit/ncenterhit));
            pwsummary.println("Average abs. activity of center bases\t"+(dsumcenterabspredict/ncentertotal));
	    //pwsummary.

	    pwsummary.close();
	    //pwstats.println("AUC\t"+szbasepredictions+"\t"+szoverlapcoordinates+"\t"+
	    //          dauc+"\t"+dauc05+"\t"+dauc01+"\t"+dauc005+"\t"+daucrev+"\t"+dauc05rev+"\t"+dauc01rev+"\t"+dauc005rev+"\t"+
	    //	daucabs+"\t"+dauc05abs+"\t"+dauc01abs+"\t"+dauc005abs+"\t"+(ntotabove/(double) (ntotabove+ntotnotabove))+"\t"+
	    //          (dsumabspredicthit/(double) ntotabove)+"\t"+(dsumabspredict/(double) (ntotabove+ntotnotabove))+"\t"+(dsumcenterabspredicthit/ncenterhit)
	    //	+"\t"+(dsumcenterabspredict/ncentertotal));
	}
	//System.out.println("SFOLD\t"+szbasepredictions+"\t"+szoverlapcoordinates+"\t"+((dsumabspredicthit/ntotabove)/(dsumabspredictnohit/ntotnotabove)));

	if (bvaluebins)
	{
	   pwstats.println("Label\tratio\tdenom\tcenterratio\ttotalratio");
	   //pwstats.println("AUC\t"+szbasepredictions+"\t"+szoverlapcoordinates+"\t"+dauc+"\t"+dauc05+"\t"+dauc01);
	   //System.out.println("SFOLD\t"+szbasepredictions+"\t"+szoverlapcoordinates+"\t"+((dsumabspredicthit/ntotabove)/(dsumabspredictnohit/ntotnotabove)));

	   int ncurrbin = 1;

           //int nlastval = baseA[0].dval;
	   double dlastval = dvaluebins*Math.round(baseA[0].dval/dvaluebins);
	   double dcurrval = 0;

	   ArrayList altotal = new ArrayList();
	   ArrayList aldenom = new ArrayList();
	   ArrayList alnumer = new ArrayList();

           for (int nindex =0; nindex < baseA.length; nindex++)
	   {
	      dcurrval = dvaluebins * Math.round(baseA[nindex].dval/dvaluebins);

	      if (dcurrval != dlastval)//||(nindex == baseA.length - 1))
	      {
	         double dratio = nsumabove/((double) nsumabove+nsumnotabove);
		 //double dcurratio = ncurrabove/((double) ncurrabove+ncurrnotabove);
		 //double dbinavg = dbinsum/nbincount;

		 altotal.add(Double.valueOf(dlastval));
		 aldenom.add(Integer.valueOf(ncurrabove+ncurrnotabove));
		 alnumer.add(Integer.valueOf(ncurrabove));

		 ncurrbin++;
		 dbinsum = 0;
		 nbincount = 0;
		 ncurrabove = 0;
		 ncurrnotabove = 0;
		 dlastval = dcurrval;
	      }

	      if (baseA[nindex].bhit)
	      {
	         nsumabove++;
		 ncurrabove++;
	      }
	      else
	      {
	         nsumnotabove++;
		 ncurrnotabove++;
	      }

	      dbinsum += baseA[nindex].dval;
	      nbincount++;
	   }

	   altotal.add(Double.valueOf(dcurrval));
	   aldenom.add(Integer.valueOf(ncurrabove+ncurrnotabove));
	   alnumer.add(Integer.valueOf(ncurrabove));

	   // int SUMTHRESH = 500;
	   int nfirstindex = 0;
	   int nsum = ((Integer) aldenom.get(nfirstindex)).intValue();

	   int nval;
           do
	   {
	      nfirstindex++;
	      nval = ((Integer) aldenom.get(nfirstindex)).intValue();
	      nsum += nval;// ((Integer) aldenom.get(nfirstindex)).intValue();
	   }
	   while ((nval < nmincountextreme)&&(nfirstindex < aldenom.size()-1));//  SUMTHRESH);

	   int nnumsum = 0;
	   int ndenomsum = 0;
	   for (int na = 0; na < nfirstindex; na++)
	   {
	      nnumsum += ((Integer) alnumer.get(na)).intValue();
	      ndenomsum += ((Integer) aldenom.get(na)).intValue();
	   }

	   int nendindex = altotal.size()-1;
	   nsum = ((Integer) aldenom.get(nendindex)).intValue();
	   do
	   {
	      nendindex--;
	      nval = ((Integer) aldenom.get(nendindex)).intValue();
              nsum += nval;// ((Integer) aldenom.get(nendindex)).intValue();
	   }
	   while ((nval < nmincountextreme)&&(nendindex>0));// SUMTHRESH);


	   int nnumsumend = 0;
           int ndenomsumend = 0;
           for (int na = nendindex+1; na < altotal.size(); na++)
	   {
	      nnumsumend += ((Integer) alnumer.get(na)).intValue();
	      ndenomsumend += ((Integer) aldenom.get(na)).intValue();
	   }

	   pwstats.println((((Double) altotal.get(nendindex+1)).doubleValue())+"\t"+(nnumsumend/(double) ndenomsumend)+"\t"+ndenomsumend+"\t"+dcenterratio+"\t"+dtotratio);
	   for (int nb = nendindex; nb >= nfirstindex; nb--)
           {
	      pwstats.println((((Double) altotal.get(nb)).doubleValue())+"\t"+(((Integer) alnumer.get(nb)).intValue())/(double) ((Integer) aldenom.get(nb)).intValue()
                            +"\t"+((Integer) aldenom.get(nb)).intValue()+"\t"+dcenterratio+"\t"+dtotratio);
	    //                            "\t"+((Integer) aldenom.get(nb)).intValue());
	   }

	   pwstats.println((((Double) altotal.get(nfirstindex-1)).doubleValue())+"\t"+(nnumsum/(double) ndenomsum)+"\t"+ndenomsum+"\t"+dcenterratio+"\t"+dtotratio);

	}
        else
	{

	    int ncurrbin = 1;

	    pwstats.println("Bin Index\tPercentile\tBin Avg\tBin Hit Fraction\tOverall Hit Fraction\tCenter Hit Fraction\t"+
                        "Cumulative Fold (vs. Overall)\tCumulative Fold Reverse Ranking (vs. Overall)");
            for (int nindex =0; nindex < baseA.length; nindex++)
	    { 
	       if (baseA[nindex].bhit)
	       {
	          nsumabove++;
	 	  ncurrabove++;
	       }
	       else
	       {
	          nsumnotabove++;
		  ncurrnotabove++;
	       }

	       if (baseA[baseA.length-nindex-1].bhit)
	       {
	          nreversesumabove++;
	       }
	       else
	       {
	          nreversesumnotabove++;
	       }

	       dbinsum += baseA[nindex].dval;
	       nbincount++;

	       if (((nindex+1)/(baseA.length/(double) numoverlapbins) >= ncurrbin))
	       {
	 	  double dratio = nsumabove/((double) nsumabove+nsumnotabove);
		  double dcurratio = ncurrabove/((double) ncurrabove+ncurrnotabove);
		  double dreverseratio = nreversesumabove/((double) nreversesumabove + nreversesumnotabove);
		  double dbinavg = dbinsum/nbincount;
		//pwstats.println((nindex+1)+"\t"+ncurrbin/(double) numoverlapbins+"\t"+dbinavg+"\t"+dratio+"\t"+dratio/dtotratio+"\t"+
		//                  dreverseratio+"\t"+dreverseratio/dtotratio+"\t"+
		//		dcurratio+"\t"+dcurratio/dtotratio);


		  pwstats.println((nindex+1)+"\t"+ncurrbin/(double) numoverlapbins+"\t"+nf.format(dbinavg)+"\t"+nf.format(dcurratio)+"\t"+nf.format(dtotratio)+"\t"+nf.format(dcenterratio)+
				"\t"+nf.format(dratio/dtotratio)+"\t"+nf.format(dreverseratio/dtotratio));


		  ncurrbin++;
	 	  dbinsum = 0;
		  nbincount = 0;
		  ncurrabove = 0;
		  ncurrnotabove = 0;
	       }
	    }
	}
	pwstats.close();
    }
 
 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static class OverlapRecCompare implements Comparator
    {

	public int compare(Object o1, Object o2)
	{
	    OverlapRec r1 = (OverlapRec) o1;
            OverlapRec r2 = (OverlapRec) o2;

	    if (r1.dval > r2.dval)
	    {
		return -1;
	    }
	    else if (r1.dval < r2.dval)
	    {
		return 1;
	    }
	    else if (r1.drandom < r2.drandom)
	    {
		return -1;
	    }
	    else if (r1.drandom > r2.drandom)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }

    public static class OverlapRec
    {
	double dval;
	boolean bhit;
	double drandom;

	OverlapRec(double dval, boolean bhit,double drandom)
	{
	    this.dval = dval;
	    this.bhit = bhit;
	    this.drandom = drandom;
	}
    }

    /**
     *
     */
    public void combinefiles() throws IOException
    {

       NumberFormat nf  = NumberFormat.getInstance();
       nf.setMaximumFractionDigits(3);

       StringTokenizer st1 = new StringTokenizer(szcombinefileset1, ",");

       PrintWriter pw = new PrintWriter(szcombineoutputfile);

	ArrayList albr1 = new ArrayList();
	while (st1.hasMoreTokens())
	{
	    albr1.add(Util.getBufferedReader(st1.nextToken()));
	}
	StringTokenizer[] stA1 = new StringTokenizer[albr1.size()];
	ArrayList albr2 = null;

	StringTokenizer[] stA2 = null;
	if (szcombinefileset2 != null)
	{
	   albr2 = new ArrayList();
	   StringTokenizer st2 = new StringTokenizer(szcombinefileset2, ",");
	   while (st2.hasMoreTokens())
           {
	       albr2.add(Util.getBufferedReader(st2.nextToken()));
	   }
	   stA2 = new StringTokenizer[albr2.size()];
	}


	String szLine;
	while ((szLine = ((BufferedReader) albr1.get(0)).readLine())!=null)
	{
	    stA1[0] = new StringTokenizer(szLine,"\t");
	    String szID = stA1[0].nextToken();
	    pw.print(szID);
	    for (int nk = 1; nk < albr1.size(); nk++)
	    {
		szLine = ((BufferedReader) albr1.get(nk)).readLine();
		if (szLine == null)
		{
                    throw new IllegalArgumentException("Not all files in "+szcombinefileset1+" have the same number of lines!");
		}
		stA1[nk] = new StringTokenizer(szLine,"\t");
		String szcurrID = stA1[nk].nextToken();
		if (!szcurrID.equals(szID))
		{
		    throw new IllegalArgumentException("We had different IDs "+szcurrID+" and "+szID+" on corresponding lines!");
		}
	    }

	    if (albr2 != null)
            {
	       for (int nk = 0; nk < albr2.size(); nk++)
	       {
		   szLine = ((BufferedReader) albr2.get(nk)).readLine();
		   if (szLine == null)
		   {
                      throw new IllegalArgumentException("Not all files in "+szcombinefileset2+" have the same number of lines with "+szcombinefileset1);
		   }
		   stA2[nk] = new StringTokenizer(szLine,"\t");
		   String szcurrID = stA2[nk].nextToken();
		   if (!szcurrID.equals(szID))
		   {
		      throw new IllegalArgumentException("We had different IDs "+szcurrID+" and "+szID+" on corresponding lines!");
		   }
	       }
	    }

	    while (stA1[0].hasMoreTokens())
	    {
		double dsum1 = Double.parseDouble(stA1[0].nextToken());
		//System.out.println("\t"+dsum1);
		for (int nk =1; nk < stA1.length; nk++)
		{
		    dsum1 += Double.parseDouble(stA1[nk].nextToken());
		    //System.out.println("\t"+dsum1);
		}

		double davg1 = dsum1/stA1.length;
		//System.out.println(davg1+"\t"+dsum1+"\t"+stA1.length);

		double doutput = davg1;

		if (albr2 != null)
		{
		   double dsum2 = Double.parseDouble(stA2[0].nextToken());

		   for (int nk =1; nk < stA2.length; nk++)
		   {
	              dsum2 += Double.parseDouble(stA2[nk].nextToken());
		   }


		   double davg2 = dsum2/stA2.length;
		   //System.out.println(doutput+"\t"+davg2);
		   if (doutput >= 0)
		   {
		       doutput = Math.max(Math.min(davg2, doutput), 0);
		   }
		   else 
		   {
		       doutput = Math.min(Math.max(davg2, doutput), 0);
		   }

		}
		pw.print("\t"+nf.format(doutput));
	    }
	    pw.println();
	}

	pw.close();
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    public void normalize() throws IOException
    {

       PrintWriter pwout = new PrintWriter(new FileWriter(sznormoutfile));
       PrintWriter pwoutnoavg = null;

       if (sznormoutfilenoavg != null)
       {
          pwoutnoavg = new PrintWriter(new FileWriter(sznormoutfilenoavg));
       }

       BufferedReader brdna = Util.getBufferedReader(szdnafile);

       String szLine;

       int numlines = 0;
       String szLine_rna;
       String szLine_dna;

       szLine_dna = brdna.readLine();
       if (szLine_dna == null)
       {
          throw new IllegalArgumentException(szdnafile+" is empty!");
       }

       StringTokenizer st = new StringTokenizer(szLine_dna,"\t");
       //the number of columns corresponds to the number of headers - assumes no header for id column
       int numcols = st.countTokens();
       double dsumdna = 0;
       double dsumrna = 0;


       //storing the probes and the total dna count
       ArrayList alprobes = new ArrayList();
       while ((szLine_dna = brdna.readLine())!=null)
       {
          st = new StringTokenizer(szLine_dna,"\t");
	  String sztoken = st.nextToken();
	  alprobes.add(sztoken);
	  while (st.hasMoreTokens())
	  {
             dsumdna += Double.parseDouble(st.nextToken())+ndnapseudocount;
	  }
	  numlines++;
       }
       brdna.close();

       brdna = Util.getBufferedReader(szdnafile);
       brdna.readLine();  //flushing the header

       BufferedReader brrna = Util.getBufferedReader(szrnafile);

       brrna.readLine();  //flushing the header

       int nline = 0;
       //double[][] dnavals = new double[numlines][numcols]; //stores the normalized DNA value
       int[][] dnacount = new int[numlines][numcols]; //stores the original count value
       int[][] rnacount = new int[numlines][numcols];

       double[][] dnaval = new double[numlines][numcols]; //stores the original count value
       double[][] rnaval = new double[numlines][numcols];
 
       double[][] ratio = new double[numlines][numcols];

       //double[][] dpredicterror = new double[NUMSTEP[MAXNEAREST+1];

       while ((szLine_rna = brrna.readLine())!=null)
       {
          szLine_dna = brdna.readLine();

	  if (szLine_dna == null)
	  {
              throw new IllegalArgumentException("RNA file "+szrnafile+" has more lines than "+szdnafile+" expecting the same number");
	  }


	  StringTokenizer strna = new StringTokenizer(szLine_rna,"\t");
	  String szrnaid = strna.nextToken();

	  StringTokenizer stdna = new StringTokenizer(szLine_dna,"\t");
	  String szdnaid = stdna.nextToken();

	  if (!szrnaid.equals(szdnaid))
	  {
	      throw new IllegalArgumentException("RNA file "+szrnafile+" does not match DNA file "+szdnafile+" with ID "+szrnaid+" "+szdnaid);
	  }
	  int ncol = 0;
	  while (strna.hasMoreTokens())
	  {
	      int nrnaval = nrnapseudocount+Integer.parseInt(strna.nextToken());
	      int ndnaval = ndnapseudocount+Integer.parseInt(stdna.nextToken());
	      dsumrna += nrnaval;
	      rnacount[nline][ncol] = nrnaval;
	      dnacount[nline][ncol] = ndnaval;
	      ncol++;
	  }
	  nline++;
       }
       brrna.close();
       brdna.close();

       double dlog2 = Math.log(2);

       for (nline = 0; nline < numlines; nline++)
       {
          for (int ncol = 0; ncol < numcols; ncol++)
	  {
	      rnaval[nline][ncol] = Math.log(rnacount[nline][ncol]/(double) dsumrna)/dlog2;

	      //System.out.println(nline+"\t"+ncol+"\t"+rnaval_rep[nline][ncol]+"\t"+rnacount_rep[nline][ncol]);
	      dnaval[nline][ncol] = Math.log(dnacount[nline][ncol]/(double) dsumdna)/dlog2;
	      ratio[nline][ncol] = rnaval[nline][ncol] - dnaval[nline][ncol];
	  }
       }

       //ArrayList[] ratios = new ArrayList[numlines];

       for (nline = 0; nline < numlines; nline++)
       {
          if (pwoutnoavg != null)
	  {
             pwoutnoavg.print(alprobes.get(nline));
	  }

	  ArrayList ratiosAL = new ArrayList();
	  for (int ncol = 0; ncol < numcols; ncol++)
	  {
	     if ((dnacount[nline][ncol]>=(ndnacutoff+ndnapseudocount)))
	     {
	        double dratioval = ratio[nline][ncol];

	        //ratios[nline].add(new Double(dratioval));
	        ratiosAL.add(Double.valueOf(dratioval));

		if (pwoutnoavg != null)
		{
		    pwoutnoavg.print("\t"+dratioval);
		}
	     }
	     else
	     {
	        if (pwoutnoavg != null)
		{
		    pwoutnoavg.print("\t");
		}
	     }
	  }

	  if (pwoutnoavg != null)
	  {
             pwoutnoavg.println();
	  }

	  //copies ratio values into array so can easily be sorted
          //double[] copyratios = new double[ratios[nline].size()];
          double[] copyratios = new double[ratiosAL.size()];
	  for (int nj = 0; nj < copyratios.length; nj++)
	  {
	      //copyratios[nj] = ((Double) ratios[nline].get(nj)).doubleValue();
             copyratios[nj] = ((Double) ratiosAL.get(nj)).doubleValue();
	  }
	  Arrays.sort(copyratios);

	  if (copyratios.length == 0)
	  {
             pwout.println(alprobes.get(nline)+"\t0\t0");
	  }
	  else if (copyratios.length % 2 == 0)
	  {
             //even number of ratio, takes the averge since in log space
	     pwout.println(alprobes.get(nline)+"\t"+(copyratios[copyratios.length/2-1]+copyratios[copyratios.length/2])/2.0+"\t"+copyratios.length);
	  }
	  else
	  {
             pwout.println(alprobes.get(nline)+"\t"+(copyratios[copyratios.length/2])+"\t"+copyratios.length);
	  }
       }

       if (pwoutnoavg != null)
       {
          pwoutnoavg.close();
       }

       pwout.close();
    }


     //////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts data into table format
     * Assumes data label, by default are of form CELL_STATE_SEQ_POSITIONINDEX_otherinfo
     * Position index starts from 0
     * Also all sequences of same row are ordered consecutively based on position being available
     */
    public void converttable() throws IOException
    {
       String szLine;
       BufferedReader br = Util.getBufferedReader(szconvertinputfile);
       PrintWriter pwconverttable = new PrintWriter(new FileWriter(szconvertoutputfile));
       int ni = 0;
       boolean bfirst = true;

       int nprevpos=-1;
       String szprevlabel="";

       HashMap hmlabels = new HashMap();
       ArrayList allabels = new ArrayList();
       int nlabelindex = 0;
       int nmaxpos = -1;
       while ((szLine = br.readLine())!=null)
       {
          StringTokenizer st = new StringTokenizer(szLine,"\t");
	  String szID = st.nextToken();
	  StringTokenizer stu = new StringTokenizer(szID,"_");

	  int npos= -1;
	  int nuindex = 0;
	  boolean bufirst = true;
	  String szlabel = "";
	  while (stu.hasMoreTokens())
	  {
	      if (nuindex == ntileindex)
	      {
	         npos = Integer.parseInt(stu.nextToken());
		 if (npos > nmaxpos)
		 {
		     nmaxpos = npos;
		 }
	      }
	      else
	      {
	         if (bufirst)
	         {
		    szlabel = stu.nextToken();
		    bufirst = false;
	         }
	         else
	         {
		     // if (nuindex <= ntileindex)
		    szlabel += "_" + stu.nextToken();
		    // else 
		    //  stu.nextToken();
		 }
	      }
	      nuindex++;
	  }
	  Integer objcount = (Integer) hmlabels.get(szlabel);
	  if (objcount == null)
	  {
	      allabels.add(szlabel);
	      hmlabels.put(szlabel, Integer.valueOf(nlabelindex));
	      nlabelindex++;
	  }
       }
       br.close();
       String[][] vals = new String[nlabelindex][nmaxpos+1];


       br = Util.getBufferedReader(szconvertinputfile);
       while ((szLine = br.readLine())!=null)
       {
          StringTokenizer st = new StringTokenizer(szLine,"\t");
	  String szID = st.nextToken();
	  StringTokenizer stu = new StringTokenizer(szID,"_");

	  int npos= -1;
	  int nuindex = 0;
	  boolean bufirst = true;
	  String szlabel = "";
	  while (stu.hasMoreTokens())
	  {
	      if (nuindex == ntileindex)
	      {
	         npos = Integer.parseInt(stu.nextToken());
	      }
	      else
	      {
	         if (bufirst)
	         {
		    szlabel = stu.nextToken();
		    bufirst = false;
	         }
	         else
	         {
		     // if (nuindex <= ntileindex)
		    szlabel += "_" + stu.nextToken();
		    // else 
		    //  stu.nextToken();
		 }
	      }
	      nuindex++;
	  }

	  int nrowindex = ((Integer) hmlabels.get(szlabel)).intValue();
          String szval = st.nextToken();
          int ncount = Integer.parseInt(st.nextToken());
	  if (ncount > 0)
	  {
	     vals[nrowindex][npos] = szval;
	  }
       }
       br.close();


       if (brecentermean|| brecenterends)
       {
	   double drecenter = 0;

	   if (brecentermean)
	   {
	       int ncount = 0;
	       double dsum = 0;

	      for (int na = 0; na < vals.length; na++)
	      {
       	         for (int nb = 0; nb < vals[na].length; nb++)
	         {
	            if (vals[na][nb] != null)
		    {
			ncount++;
			dsum += Double.parseDouble(vals[na][nb]);
		    }
		 }
	      }

	      drecenter = dsum/ncount;
	   }
	   else
	   {
               int ncount = 0;
               double dsum = 0;

	       for (int na = 0; na < vals.length; na++)
	       {
	          if (vals[na][0] != null)
	       	  {
		      ncount++;
		      dsum += Double.parseDouble(vals[na][0]);
		  }	

                  if (vals[na][vals[na].length-1] != null)
	          {
		      ncount++;
		      dsum += Double.parseDouble(vals[na][vals[na].length-1]);
		  }
       	   
	       }
	       drecenter = dsum/ncount;

	   }


	   for (int na = 0; na < vals.length; na++)
	   {
       	      pwconverttable.print(allabels.get(na));
       	      for (int nb = 0; nb < vals[na].length; nb++)
	      {
	         if (vals[na][nb] == null)
		 {
		     pwconverttable.print("\t");
		 }
	         else
	         {
		     pwconverttable.print("\t"+(Double.parseDouble(vals[na][nb])-drecenter));
		 }
	      }
       	      pwconverttable.println();
	   }
	   pwconverttable.close();
       }
       else
       {

          for (int na = 0; na < vals.length; na++)
          {
	     pwconverttable.print(allabels.get(na));
	     for (int nb = 0; nb < vals[na].length; nb++)
	     {
	        if (vals[na][nb] == null)
	        {
		   pwconverttable.print("\t");
	        }
	        else
	        {
                   pwconverttable.print("\t"+vals[na][nb]);
	        }
	     }
	     pwconverttable.println();
	  }
          pwconverttable.close();
       }

    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Converts data into table format
     * Assumes data label, by default are of form CELL_STATE_SEQ_POSITIONINDEX_otherinfo
     * Position index starts from 0
     * Also all sequences of same row are ordered consecutively based on position being available
     */
/*
    public void converttable() throws IOException
    {
       String szLine;
       BufferedReader br = new BufferedReader(new FileReader(szconvertinputfile));
       PrintWriter pwconverttable = new PrintWriter(new FileWriter(szconvertoutputfile));
       int ni = 0;
       boolean bfirst = true;

       int nprevpos=-1;
       String szprevlabel="";

       while ((szLine = br.readLine())!=null)
       {
          StringTokenizer st = new StringTokenizer(szLine,"\t");
	  String szID = st.nextToken();
	  StringTokenizer stu = new StringTokenizer(szID,"_");

	  int npos= -1;
	  int nuindex = 0;
	  boolean bufirst = true;
	  String szlabel = "";
	  while (stu.hasMoreTokens())
	  {
	      if (nuindex == ntileindex)
	      {
	         npos = Integer.parseInt(stu.nextToken());
	      }
	      else
	      {
	         if (bufirst)
	         {
		    szlabel = stu.nextToken();
		    bufirst = false;
	         }
	         else
	         {
		     // if (nuindex <= ntileindex)
		    szlabel += "_" + stu.nextToken();
		    // else 
		    //  stu.nextToken();
		 }
	      }
	      nuindex++;
	  }
	
	  if (!szlabel.equals(szprevlabel))
	  {
	      //we are on a new row
	     if (!bfirst)
	     {
		 //if not first row ends previous line
	        pwconverttable.println();
	     }

	     pwconverttable.print(szlabel);

	     //previous position is before index 0
	     nprevpos = -1;
	  }

	  for (int nj = nprevpos+1; nj < npos; nj++)
	  {
	      //fills in missing tabs up to current position
	     pwconverttable.print("\t");
	  }

	  //gets the actual value associated with this position
	  String szval = st.nextToken();
	  int ncount = Integer.parseInt(st.nextToken());
	  if (ncount == 0)
	  {
	      //szval not based on any measures so skipping
	      pwconverttable.print("\t");
	  }
          else
	  {
	      //reporting szval value
	      pwconverttable.print("\t"+szval);
	  }

	  //update the previous values
	  bfirst = false;
	  szprevlabel = szlabel;
	  nprevpos = npos;
       }

       pwconverttable.println();
       pwconverttable.close();
       br.close();
    }
*/

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    /**
     * step should be less than tile length
     */ 
    public void deconvolveExact() throws IOException
    {
       NumberFormat nf  = NumberFormat.getInstance();
       nf.setMaximumFractionDigits(3);

       if (ntilelength % nstepsize != 0)
	   throw new IllegalArgumentException("tilelength is not divisible by stepsize values are "+ntilelength+" "+nstepsize);
 
       int numhiddenvalstotal = ((ntilelength) + nstepsize * (numtiles-1))/nstepsize;



       //System.out.println("num hiddenvals = "+numhiddenvalstotal+" ngcd = "+ngcd);

       //int NUMLINES = 7860;
       //int NUMCOLS = 31;
       //int NUMHIDDENVALS = 59;
       //int NUMHIDDENCOLS = 5;

       //int nseqlength = ntilelength/ngcd;
       int numhiddenvalstile = ntilelength/nstepsize;
       //int SEQLENGTH = 29;
       //int numhiddenvalstilesq = numhiddenvalstile * numhiddenvalstile;
       //int SEQLENGTHSQ = SEQLENGTH * SEQLENGTH;

       //Random theRandom = new Random(4343);


       String szLine;

       BufferedReader brtable;

       brtable =  Util.getBufferedReader(szdeconvolveinputfile);
       int numlines = 0;

       while ((szLine = brtable.readLine())!=null)
       {
	   numlines++;
       }
       brtable.close();
        

       brtable =  Util.getBufferedReader(szdeconvolveinputfile);
       PrintWriter pw = new PrintWriter(new FileWriter(szdeconvolveoutputfile));

       double[][] data = new double[numlines][numtiles];
       boolean[][] bpresent = new boolean[numlines][numtiles];
       String[] labels = new String[numlines];

       //ArrayList alvals = new ArrayList();
       int nrow = 0;
       int npos = 0;
       double dval = 0;
       double dvalsq = 0;
 
       while ((szLine = brtable.readLine())!=null)
       {
          StringTokenizer st = new StringTokenizer(szLine,"\t",true);
	  String szid = st.nextToken();
	  st.nextToken();
	  labels[nrow] = szid;
	  int ncol = 0;
	  while (st.hasMoreTokens())
	  {
	     if (ncol >= numtiles)
	     {
		 throw new IllegalArgumentException("We have more entries in "+szdeconvolveinputfile+" for "+szid+" than number of tile positions specified "+numtiles);
	     }

	     String sztoken = st.nextToken();
	     if (!sztoken.equals("\t"))
	     {
		 data[nrow][ncol] = Double.parseDouble(sztoken);
		 bpresent[nrow][ncol] = true;
		 if (st.hasMoreTokens())
		 {
		    st.nextToken();
		 }

		 npos++;
		 dval += data[nrow][ncol];
		 dvalsq += data[nrow][ncol]*data[nrow][ncol];
	     }
	     else
	     {
		 bpresent[nrow][ncol] = false;
		 data[nrow][ncol] =0;
	     }
	     ncol++;
	  }

	  while (ncol < numtiles)
	  {
	      bpresent[nrow][ncol] = false;
	      data[nrow][ncol] =0;
	      ncol++;
	  }
	  nrow++;
       }
       brtable.close();
   
       double dinitmean = dval/npos;
       double dnormvar = dvalsq/npos -(dval/npos)*(dval/npos);

       /*
       boolean BZSCORE = false;//true;

       if (BZSCORE)
       {
	   double dstd = Math.sqrt(dnormvar);

	   for (nrow = 0; nrow < data.length; nrow++)
	   {
	       for (int ncol = 0; ncol < data[nrow].length; ncol++)
	       {
		   if (bpresent[nrow][ncol])
		   {
		       data[nrow][ncol] = (data[nrow][ncol]-dinitmean)/dstd;
		   }
	       }
	   }
	   System.out.println("normalizing z score\t"+dinitmean+"\t"+dstd);
	   dinitmean = 0;
	   dnormvar = 1;
       }
       */

       //double[] datavals = new double[numhiddenvalstotal];

       //double[] datavalssum = new double[numhiddenvalstotal];
       double[][] dhiddenvalA = new double[numlines][numhiddenvalstotal];

       double[][] meanX = new double[numhiddenvalstotal][1];

       for (int nj = 0; nj < numhiddenvalstotal; nj++)
       {
          meanX[nj][0] = dinitmean;
       }


       RealMatrix theRealMatrixMeanX = MatrixUtils.createRealMatrix(meanX);


       for (nrow = 0; nrow < data.length; nrow++)
       {
	   //System.out.println(nrow);
	   //for (int nval = 0; nval < datavals.length; nval++)
	   //{
	   //  datavals[nval] = dinitmean;
	   //  datavalssum[nval] = 0;
	   //}
	  double[] data_nrow = data[nrow];
	  boolean[] bpresent_nrow = bpresent[nrow];

	  int numpresenttiles = 0;
	  for (int nk = 0; nk < bpresent_nrow.length; nk++)
	  {
	      if (bpresent_nrow[nk])
	      {
		  numpresenttiles++;
	      }
	  }

	  if (numpresenttiles == 0)
	  {
	      for (int nk = 0; nk < dhiddenvalA[nrow].length; nk++)
	      {
	         dhiddenvalA[nrow][nk] = dinitmean;
	      }
	  }
	  else
	  {

             double[][] covarXY = new double[numhiddenvalstotal][numpresenttiles];
             double[][] covarYY = new double[numpresenttiles][numpresenttiles]; 
             double[][] observedY = new double[numpresenttiles][1];
             double[][] meanY = new double[numpresenttiles][1];

	     int[] presentindiciesA = new int[numpresenttiles];

	     int npresenttile = 0;
	     for (int ntile = 0; ntile < numtiles; ntile++)
	     {
	        if (bpresent_nrow[ntile])
	        {
		   observedY[npresenttile][0] = data_nrow[ntile];
		   meanY[npresenttile][0] = dinitmean;

		   presentindiciesA[npresenttile] = ntile;


		   for (int nxval = ntile; nxval < ntile + numhiddenvalstile; nxval++)
		   {
	              covarXY[nxval][npresenttile] = dpriorvar/(double) numhiddenvalstile;
		   } 

		   npresenttile++;
		}
	     }


	     for (npresenttile = 0; npresenttile < numpresenttiles; npresenttile++)
	     {
	        for (int npresenttile2 = 0; npresenttile2 < numpresenttiles; npresenttile2++)
	        {
		   if (npresenttile == npresenttile2)
		   {
		       covarYY[npresenttile][npresenttile2] = dnormvar + dpriorvar/(double) numhiddenvalstile;
		   }
		   else
		   {
		      if (Math.abs(presentindiciesA[npresenttile]-presentindiciesA [npresenttile2]) < numhiddenvalstile)
		      {
			  covarYY[npresenttile][npresenttile2] = dpriorvar * (numhiddenvalstile-Math.abs(presentindiciesA[npresenttile]-presentindiciesA[npresenttile2]))/
                                                                                     (double) (numhiddenvalstile*numhiddenvalstile);
		      } 
		   }
		}
	     }
	 

	     RealMatrix theRealMatrixMeanY = MatrixUtils.createRealMatrix(meanY);
             RealMatrix theRealMatrixObservedY = MatrixUtils.createRealMatrix(observedY);
             RealMatrix theRealMatrixCovarXY = MatrixUtils.createRealMatrix(covarXY);
             RealMatrix theRealMatrixCovarYY = MatrixUtils.createRealMatrix(covarYY);

             dhiddenvalA[nrow] = theRealMatrixMeanX.add(theRealMatrixCovarXY.multiply(MatrixUtils.inverse(theRealMatrixCovarYY).multiply(							                                                    theRealMatrixObservedY.subtract(theRealMatrixMeanY)))).getColumn(0);
	     /*
	     if (nrow == 0)
	     {
		 System.out.println("meanY");
		 for (int na = 0; na < meanY.length; na++)
		 {
		     for (int nb = 0; nb < meanY[na].length; nb++)
		     {
			 System.out.print(meanY[na][nb]+"\t");
		     }
		     System.out.println();
		 }


		 System.out.println("observedY");
		 for (int na = 0; na < observedY.length; na++)
		 {
		     for (int nb = 0; nb < observedY[na].length; nb++)
		     {
			 System.out.print(observedY[na][nb]+"\t");
		     }
		     System.out.println();
		 }


	         System.out.println("meanX");
		 for (int na = 0; na < meanX.length; na++)
		 {
		     for (int nb = 0; nb < meanX[na].length; nb++)
		     {
			 System.out.print(meanX[na][nb]+"\t");
		     }
		     System.out.println();
		 }


	         System.out.println("covarXY");
		 for (int na = 0; na < covarXY.length; na++)
		 {
		     for (int nb = 0; nb < covarXY[na].length; nb++)
		     {
			 System.out.print(covarXY[na][nb]+"\t");
		     }
		     System.out.println();
		 }



	         System.out.println("covarYY");
		 for (int na = 0; na < covarYY.length; na++)
		 {
		     for (int nb = 0; nb < covarYY[na].length; nb++)
		     {
			 System.out.print(covarYY[na][nb]+"\t");
		     }
		     System.out.println();
		 }

	     }
	     */
	
	  }
       }


       if (bstandardize)
       {
          double dsumhidden = 0;
          double dsumhiddensq = 0;
	  int numels = 0;
	  for (int nrowindex = 0; nrowindex < dhiddenvalA.length; nrowindex++)
	  {
	      for (int nk = 0; nk < dhiddenvalA[nrowindex].length; nk++)
	      {
		  double dcurrval = dhiddenvalA[nrowindex][nk];
		  dsumhidden += dcurrval;
		  dsumhiddensq += dcurrval * dcurrval;
		  numels++;
	      }
	  }

	  double dmeanhidden = dsumhidden/numels;
	  double dvarhidden = dsumhiddensq/numels - dmeanhidden * dmeanhidden;
	  double dstdhidden = Math.sqrt(dvarhidden);

          for (int nrowindex = 0; nrowindex < labels.length; nrowindex++)
          {
	     pw.print(labels[nrowindex]);
	     double[] dhiddenvalA_nrowindex = dhiddenvalA[nrowindex];
	     for (int nk = 0; nk < dhiddenvalA_nrowindex.length; nk++)
             {
		 pw.print("\t"+nf.format((dhiddenvalA_nrowindex[nk]-dmeanhidden)/dstdhidden));
	     }
	     pw.println();
	  }
	  pw.close();
       }
       else
       {
          for (int nrowindex = 0; nrowindex < labels.length; nrowindex++)
          {
	     pw.print(labels[nrowindex]);
             double[] dhiddenvalA_nrowindex = dhiddenvalA[nrowindex];
	     for (int nk = 0; nk < dhiddenvalA_nrowindex.length; nk++)
             {
	        pw.print("\t"+nf.format(dhiddenvalA_nrowindex[nk]));
	     }
	     pw.println();
	  }
	  pw.close();
       }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void linearInterpolate() throws IOException
    {
	// this.szinterpolateinputfile = szinterpolateinputfile;
        //this.szinterpolateoutputfile = szinterpolateoutputfile;

	NumberFormat nf  = NumberFormat.getInstance();
	nf.setMaximumFractionDigits(3);

	BufferedReader br = Util.getBufferedReader(szinterpolateinputfile);
	PrintWriter pw = new PrintWriter(new FileWriter(szinterpolateoutputfile));
	String szLine;

	double[] data = null;
	while ((szLine = br.readLine())!=null)
	{
	   StringTokenizer st = new StringTokenizer(szLine,"\t");
	   if (data == null)
	   {
	      data = new double[st.countTokens()-1];
	   }

	   pw.print(st.nextToken());

	   int nindex = 0;
	   while (st.hasMoreTokens())
           {
	      data[nindex] = Double.parseDouble(st.nextToken());
	      nindex++;
	   }

	   int numpos = nstepsize*data.length;
  	   for (int npos = 0; npos < numpos; npos++)
	   {
	      double dsmoothval;

	      int nhiddenindex = npos/nstepsize;
	      int nhiddenpos = npos % nstepsize;

	      //4

              if (((nhiddenindex ==0)&&(nhiddenpos <nstepsize/2))||(nhiddenpos == nstepsize/2)||((nhiddenindex == data.length-1)&&(nhiddenpos>=nstepsize/2)))
	      {
	         dsmoothval = data[nhiddenindex];
	      }
	      else if (nhiddenpos < nstepsize/2)
	      {
		  //0,1
		  //0--> 3/4, 1/4
		  //1-->1, 0
	         dsmoothval = ((nstepsize/2+nhiddenpos+1)/(double) nstepsize)*data[nhiddenindex]+(1-(nstepsize/2+nhiddenpos+1)/(double) nstepsize)*data[nhiddenindex-1];
	      }
	      else
	      {
		  //2-->1,0
		  //3-->3/4,1/4
	         dsmoothval = (1-(nhiddenpos-nstepsize/2)/(double) nstepsize)*data[nhiddenindex]+((nhiddenpos-nstepsize/2)/(double) nstepsize)*data[nhiddenindex+1];
	      }

	      pw.print("\t"+nf.format(dsmoothval));
	   }
	   pw.println();
	}
	br.close();
	pw.close();
    }


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public void adjacentchanges() throws IOException
    {
	//param for num entries per line
	//double dpval;
	StringTokenizer stdatafiles = new StringTokenizer(szdatafiles,",");
	int nfile = 0;

	PrintWriter pw = new PrintWriter(new FileWriter(szadjacentoutfile));

	NumberFormat nf = NumberFormat.getInstance();
	nf.setMaximumFractionDigits(6);

	int numfiles = stdatafiles.countTokens();
	ArrayList[][][] alvals = null;
	String[][][] ids = null;
	double[][][] pvalsLESS = null;
        double[][][] pvalsGREATER = null;
	//double[][][] pvals = null;

	//double[][][] randomizedpvalsLESS = null;
	//double[][][] randomizedpvalsGREATER = null;



	double[][] combinedpvalsLESS = null;
        double[][] combinedpvalsGREATER = null;

	//double[][] combinedpvals = null;

        //double[][] combinedrandomizedpvalsLESS = null;
        //double[][] combinedrandomizedpvalsGREATER = null;

	//double[][] combinedrandomizedpvals = null;
	double[][][] medians = null;

	//double[][][] randomizedmedians = null;

	String szLine;



	while (stdatafiles.hasMoreTokens())
	{
	    String szdatafile = stdatafiles.nextToken();

	    if (nfile == 0)
	    {
	       BufferedReader br = Util.getBufferedReader(szdatafile);

	       int numlines = 0;
	       while ((szLine = br.readLine())!=null)
	       {
	          numlines++;
	       }	    
	       br.close();

	       if (numlines % numtiles != 0)
	       {
	          throw new IllegalArgumentException("The total number of lines in "+szdatafile+" is not divisible by the number of positions per regulatory element "+numtiles);
	       }

	       int numseqs = numlines/numtiles;

	       alvals = new ArrayList[numfiles][numseqs][numtiles];
	       ids = new String[numfiles][numseqs][numtiles];
	       pvalsLESS = new double[numfiles][numseqs][numtiles];
	       pvalsGREATER = new double[numfiles][numseqs][numtiles];
	       //pvals = new double[numfiles][numseqs][numtiles];
	       //randomizedpvalsLESS = new double[numfiles][numseqs][numtiles];
	       //randomizedpvalsGREATER = new double[numfiles][numseqs][numtiles];

	       combinedpvalsLESS = new double[numseqs][numtiles];
               combinedpvalsGREATER = new double[numseqs][numtiles];

               //combinedpvals = new double[numseqs][numtiles];

	       //combinedrandomizedpvalsLESS = new double[numseqs][numtiles];
	       //combinedrandomizedpvalsGREATER = new double[numseqs][numtiles];

	       medians = new double[numfiles][numseqs][numtiles];
	       //randomizedmedians = new double[numfiles][numseqs][numtiles];

	    }

	    BufferedReader br = Util.getBufferedReader(szdatafile);

	    int nseq = 0;
	    while ((szLine = br.readLine())!=null)
	    {
		//System.out.println(szLine);
		//System.out.println(nseq+"\t"+alvals.length+"\t"+alvals[nfile].length+"\t");
		for (int ntile = 0; ntile < alvals[nfile][nseq].length; ntile++)
		{
		    if (ntile >=1)
		    {
		       szLine = br.readLine();
		       if (szLine == null)
		       {
			   throw new IllegalArgumentException("In "+szdatafile+" we don't have data for every sequence and tile position expecting");
		       }
		    }
		    alvals[nfile][nseq][ntile] = new ArrayList();
		    StringTokenizer st =new StringTokenizer(szLine,"\t");
		    //st.nextToken();
		    ids[nfile][nseq][ntile] = st.nextToken();
		    if (!ids[0][nseq][ntile].equals(ids[nfile][nseq][ntile]))
		    {
			throw new IllegalArgumentException("We have mismatching IDs in the data files. "+ids[0][nseq][ntile]+" does not match "+ids[nfile][nseq][ntile]);
		    }

		    while (st.hasMoreTokens())
		    {
	               alvals[nfile][nseq][ntile].add(Double.valueOf(Double.parseDouble(st.nextToken())));		       
		    }
		}
	       
	       nseq++;
	    }
	    br.close();
	    nfile++;
	}





	//ArrayList[][][] randomizedvals = new ArrayList[alvals.length][alvals[0].length][alvals[0][0].length];
	//Random theRandom = new Random(563);
	//createrandomized(alvals,randomizedvals, theRandom);

	computepvalues(alvals, pvalsLESS,pvalsGREATER,combinedpvalsLESS,combinedpvalsGREATER, medians);

	//computepvalues(alvals, pvalsLESS, combinedpvalsLESS,pvalsGREATER,combinedpvalsGREATER,medians);
	//computepvalues(randomizedvals, randomizedpvalsLESS, combinedrandomizedpvalsLESS,randomizedpvalsGREATER, combinedrandomizedpvalsGREATER, randomizedmedians);

        double[][] combinedpvalsmerged = new double[combinedpvalsLESS.length][combinedpvalsLESS[0].length];
        double[] combinedpvalsALL = new double[combinedpvalsLESS.length*(combinedpvalsLESS[0].length-1)];
	double[] fdrs = new double[combinedpvalsALL.length];

        int nd = 0;
        for (int na = 0; na < combinedpvalsLESS.length; na++)
        {
	   for (int nb = 1; nb < combinedpvalsLESS[na].length; nb++)
           {
	       double dpval= Math.min(1,2*Math.min(combinedpvalsLESS[na][nb],combinedpvalsGREATER[na][nb]));
	       combinedpvalsmerged[na][nb] = dpval;
	       combinedpvalsALL[nd] = dpval;
	       nd++;
	   }
	}


	/*
	double[] fullrandomizedpvals = new double[2*randomizedpvalsLESS[0].length*(randomizedpvalsLESS[0][0].length-1)];

	int nd = 0;
	for (int na = 0; na < combinedrandomizedpvalsLESS.length; na++)
	{
           for (int nb = 1; nb < combinedrandomizedpvalsLESS[na].length; nb++)
	   {
       	      fullrandomizedpvals[nd] = combinedrandomizedpvalsLESS[na][nb];
	      nd++;
           }
	}

        for (int na = 0; na < combinedrandomizedpvalsGREATER.length; na++)
	{
       	   for (int nb = 1; nb < combinedrandomizedpvalsGREATER[na].length; nb++)
           {
	       fullrandomizedpvals[nd] = combinedrandomizedpvalsGREATER[na][nb];
	       nd++;
	   }
	}


	Arrays.sort(fullrandomizedpvals);

	double[] correctedpvalsALL = new double[fullrandomizedpvals.length];
	double[][] correctedpvalsGREATER = new double[combinedpvalsGREATER.length][combinedpvalsGREATER[0].length];
	double[][] correctedpvalsLESS = new double[combinedpvalsLESS.length][combinedpvalsLESS[0].length];
	double[] fdrs = new double[fullrandomizedpvals.length];

	nd = 0;
	for (int na = 0; na < combinedpvalsGREATER.length; na++)
	{
	   for (int nb = 1; nb < combinedpvalsGREATER[na].length; nb++)
       	   {
	       int nrandindex= Arrays.binarySearch(fullrandomizedpvals, combinedpvalsGREATER[na][nb]);
	       if (nrandindex < 0)
	       {
	          nrandindex = -nrandindex - 1;
	       }
	       correctedpvalsGREATER[na][nb] = (nrandindex+1)/(double) fullrandomizedpvals.length;
	       correctedpvalsALL[nd] = (nrandindex+1)/(double) fullrandomizedpvals.length;
	       nd++;
	   }
	}


        for (int na = 0; na < combinedpvalsLESS.length; na++)
	{
	   for (int nb = 1; nb < combinedpvalsLESS[na].length; nb++)
           {
	      int nrandindex= Arrays.binarySearch(fullrandomizedpvals, combinedpvalsLESS[na][nb]);
	      if (nrandindex < 0)
	      {
	         nrandindex = -nrandindex - 1;
	      }

	      correctedpvalsLESS[na][nb] = (nrandindex+1)/(double) fullrandomizedpvals.length;
	      correctedpvalsALL[nd] = (nrandindex+1)/(double)fullrandomizedpvals.length;
	      nd++;
	   }
	}

	*/

	Arrays.sort(combinedpvalsALL);

	for (int na = 0; na < fdrs.length; na++)
	{
	    fdrs[na] = combinedpvalsALL[na]*combinedpvalsALL.length/(na+1);
	    //System.out.println(na+"\t"+correctedpvals[na]+"\t"+fdrs[na]);
	}

	for (int na = fdrs.length-2; na >= 0; na--)
	{
	    fdrs[na] = Math.min(fdrs[na],fdrs[na+1]);
	}
	

	//for (int ni = 0; ni < fullrandomizedpvals.length; ni++)
	//  System.out.println(ni+"\t"+fullrandomizedpvals[ni]);

	//this.szsequences = szsequences;
	//this.nextension = nextension;
	BufferedReader brsequences = Util.getBufferedReader(szsequences);
	HashMap hmsequences = new HashMap();

	while ((szLine = brsequences.readLine())!=null)
	{
	    StringTokenizer st = new StringTokenizer(szLine,"\t");
	    hmsequences.put(st.nextToken(),st.nextToken());
	}
	brsequences.close();

	HashMap hmcoords = new HashMap();

        BufferedReader brcoords = Util.getBufferedReader(szreportercoordinates);
        while ((szLine = brcoords.readLine())!=null)
        {
	   StringTokenizer st = new StringTokenizer(szLine,"\t");
      	   hmcoords.put(st.nextToken(),st.nextToken()+"\t"+st.nextToken()+"\t"+st.nextToken());
        }
        brcoords.close();


	pw.print("ReporterID1\tReporterID2\t");
	pw.print("\tp-val_(reporter1<>reporter2)\tfdr_(reporter1<>reporter2)\tactivating_reporterID\trepressive_reporterID");
	pw.print("\tactivating_seq\trepressive_seq\tactivating_chrom\tactivating_start\tactivating_end\trepressive_chrom\trepressive_start\trepressive_end");
	for (int nc = 0; nc < pvalsLESS.length; nc++)
	{
	    pw.print("\tp-val_(reporter2<reporter1)_exp"+(nc+1));
	}

        for (int nc = 0; nc < pvalsGREATER.length; nc++)
	{
	    pw.print("\tp-val_(reporter1<reporter2)_exp"+(nc+1));
        }
	pw.print("\tp-val_(reporter2<reporter1)_ALL\tp-val_(reporter1<reporter2)_ALL");
        for (int nc = 0; nc < medians.length; nc++)
	{
	    pw.print("\tmedian_reporter1_exp"+(nc+1));
        }

        for (int nc = 0; nc < medians.length; nc++)
	{
	    pw.print("\tmedians_reporter2_exp"+(nc+1));
        }


        pw.println();




        for (int na = 0; na < pvalsGREATER[0].length; na++)
	{
	   for (int nb = 1; nb < pvalsGREATER[0][na].length; nb++)
           {
	       pw.print(ids[0][na][nb-1]+"\t"+ids[0][na][nb]);//na+"\t"+nb);

	       //double dcorrectedpval = combinedpvals[na][nb] * fullrandomizedpvals.length;
	       //double dcorrectedpvalGREATER = combinedpvalsGREATER[na][nb] * 2*pvalsGREATER[0].length*(pvalsGREATER[0][0].length-1);
	       //double dcorrectedpvalLESS = combinedpvalsLESS[na][nb] * 2*pvalsLESS[0].length*(pvalsLESS[0][0].length-1);
		   //int nrandindex= Arrays.binarySearch(fullrandomizedpvals, combinedpvals[na][nb]);
		   //if (nrandindex < 0)
		   //nrandindex = -nrandindex - 1;

		   //double dfdr = nrandindex/(double) fullrandomizedpvals.length;
	       int ncombinedpvalindex = Arrays.binarySearch(combinedpvalsALL, combinedpvalsmerged[na][nb]);
	       String szseq1 = (String) hmsequences.get(ids[0][na][nb-1]);
	       String szseq2 = (String) hmsequences.get(ids[0][na][nb]);
	       String szcoord1 = (String) hmcoords.get(ids[0][na][nb-1]);
	       String szcoord2 = (String) hmcoords.get(ids[0][na][nb]);

	       StringTokenizer st1 = new StringTokenizer(szcoord1,"\t");
	       StringTokenizer st2 = new StringTokenizer(szcoord2,"\t");

	       String szchrom1 = st1.nextToken();
	       int nbegin1 = Integer.parseInt(st1.nextToken());
	       //int nend1 =
               Integer.parseInt(st1.nextToken());


	       // System.out.println("****"+szcoord2);
	       String szchrom2 = st2.nextToken();
	       //int nbegin2 = 
               Integer.parseInt(st2.nextToken());
	       int nend2 = Integer.parseInt(st2.nextToken());


	       pw.print("\t"+combinedpvalsmerged[na][nb]+"\t"+fdrs[ncombinedpvalindex]);

	       if (combinedpvalsLESS[na][nb] <= combinedpvalsGREATER[na][nb])
	       {
	          pw.print("\t"+ids[0][na][nb-1]+"\t"+ids[0][na][nb]+"\t"+
                                      szseq1.substring(0,nstepsize+nextension)+"\t"+szseq2.substring(szseq2.length()-(nstepsize+nextension),szseq2.length())+
				     "\t"+szchrom1+"\t"+nbegin1+"\t"+(nbegin1+nstepsize+nextension)+"\t"+szchrom2+"\t"+(nend2-(nstepsize+nextension))+"\t"+nend2);
	       }
	       else
	       {
		  pw.print("\t"+ids[0][na][nb]+"\t"+ids[0][na][nb-1]+"\t"+
                               szseq2.substring(szseq2.length()-(nstepsize+nextension),szseq2.length())+"\t"+szseq1.substring(0,nstepsize+nextension)+
				      "\t"+szchrom2+"\t"+(nend2-(nstepsize+nextension))+"\t"+nend2+"\t"+szchrom1+"\t"+nbegin1+"\t"+
				      (nbegin1+nstepsize+nextension));
	       }

               for (int nc = 0; nc < pvalsLESS.length; nc++)
	       {
	          pw.print("\t"+pvalsLESS[nc][na][nb]);
	       }
	       for (int nc = 0; nc < pvalsGREATER.length; nc++)
	       { 
		  pw.print("\t"+pvalsGREATER[nc][na][nb]);
	       }

	       pw.print("\t"+combinedpvalsLESS[na][nb]+"\t"+combinedpvalsGREATER[na][nb]);               
               for (int nc = 0; nc < medians.length; nc++)
	       {
	          pw.print("\t"+nf.format(medians[nc][na][nb-1]));
	       }
               for (int nc = 0; nc < medians.length; nc++)
 	       {
	          pw.print("\t"+nf.format(medians[nc][na][nb]));
	       }
	       pw.println();



	       
	   }
	}


	pw.close();

	/*

        for (int na = 0; na < pvalsLESS[0].length; na++)
	{
	   for (int nb = 1; nb < pvalsLESS[0][na].length; nb++)
           {
	       System.out.print(na+"\t"+nb);
	       for (int nc = 0; nc < pvalsLESS.length; nc++)
	       { 
		   System.out.print("\t"+pvalsLESS[nc][na][nb]);
	       }
               for (int nc = 0; nc < pvalsLESS.length; nc++)
	       {
	           System.out.print("\t"+nf.format(medians[nc][na][nb]));
	       }
               for (int nc = 0; nc < pvalsLESS.length; nc++)
	       {
	           System.out.print("\t"+nf.format(medians[nc][na][nb-1]));
	       }


	       double dcorrectedpvalLESS = combinedpvalsLESS[na][nb] * 2*pvalsLESS[0].length*(pvalsLESS[0][0].length-1);
		   //int nrandindex= Arrays.binarySearch(fullrandomizedpvals, combinedpvals[na][nb]);
		   //if (nrandindex < 0)
		   //nrandindex = -nrandindex - 1;

		   //double dfdr = nrandindex/(double) fullrandomizedpvals.length;
	       String szseq1 = (String) hmsequences.get(ids[0][na][nb-1]);
	       String szseq2 = (String) hmsequences.get(ids[0][na][nb]);
	       int ncorrectedpvalindex = Arrays.binarySearch(combinedpvalsALL, combinedpvalsLESS[na][nb]);

	       System.out.print("\t"+combinedpvals[na][nb]+"\t"+Math.min(dcorrectedpvalLESS,1)+"\t"+correctedpvalsLESS[na][nb]+"\t"+fdrs[ncorrectedpvalindex]);

	       //System.out.println("\t"+szseq2.substring(szseq2.length()-(nstepsize+nextension),szseq2.length())+"\t"+szseq1.substring(0,nstepsize+nextension));
	   }
	}
	*/



	/*
	double dtest = -2*(Math.log(data[ncol])+
			   Math.log(data2[ncol]));
	ChiSquaredDistribution theChiSq = new ChiSquaredDistribution(4);
	double dpval =  1-theChiSq.cumulativeProbability(dtest);
	pvalA[nindex] = 1-dpval;
	nindex++;
	pvalA[nindex] = -(1-dpval);
	nindex++;
	*/

	// dpval = theMannWhitneyUTest.mannWhitneyUTest(copyratios,prevcopyratios);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    private void computepvalues(ArrayList[][][] alvals, double[][][] pvalsLESS, double[][][] pvalsGREATER, double[][] combinedpvalsLESS,double[][] combinedpvalsGREATER, 
				double[][][] medians)
    {
	
	for (int nfile = 0; nfile < alvals.length; nfile++)
	{
	    for (int nseq = 0; nseq < alvals[nfile].length; nseq++)
	    {
		double[] prevratios = new double[alvals[nfile][nseq][0].size()];
		for (int nk = 0; nk < prevratios.length; nk++)
		{
		    prevratios[nk] = ((Double) alvals[nfile][nseq][0].get(nk)).doubleValue();
		}
		Arrays.sort(prevratios);

		if (prevratios.length ==0)
		{
		    medians[nfile][nseq][0] = 0;
		}
		else if (prevratios.length % 2 == 1)
                {
                   medians[nfile][nseq][0] = prevratios[prevratios.length/2];
                }
		else
                {
                   medians[nfile][nseq][0] = (prevratios[prevratios.length/2]+prevratios[prevratios.length/2-1])/2;
                }


		for (int npos = 1; npos < alvals[nfile][nseq].length; npos++)
		{
		    double[] curratios = new double[alvals[nfile][nseq][npos].size()];
		    for (int nk = 0; nk < curratios.length; nk++)
		    {
		       curratios[nk] = ((Double) alvals[nfile][nseq][npos].get(nk)).doubleValue();
		    }
		    MannWhitneyUTest theMannWhitneyUTest = new MannWhitneyUTest();
		    double dpval = 1;

		    RankRec[] rankRecA = new RankRec[curratios.length+prevratios.length];
		    int nb = 0;
		    for (int na = 0; na < curratios.length; na++)
		    {
			rankRecA[nb] = new RankRec(curratios[na], false);
			nb++;
		    }


                    for (int na = 0; na < prevratios.length; na++)
		    {
		        rankRecA[nb] = new RankRec(prevratios[na], true);
		        nb++;
		    }

		    Arrays.sort(rankRecA, new RankRecCompare());

		    int nlastrank = 0;
		    double davgcurr = 0;
		    double davgprev = 0;

		    for (int nk = 0; nk < rankRecA.length; nk++)
		    {
			if ((nk+1 == rankRecA.length)||(rankRecA[nk].dval != rankRecA[nk+1].dval))
			{
			    int nassignrank = (nk+nlastrank)/2;
			    for (int nj = nlastrank; nj <= nk; nj++)
			    {
				if (rankRecA[nk].bprev)
			        {
				    davgprev += nassignrank;//rankRecA[nk].dval;
				}
				else
				{
				    davgcurr += nassignrank;//rankRecA[nk].dval;
				}
			    }
			    nlastrank = nk+1;
			}
		    }

		    davgprev /= prevratios.length;
		    davgcurr /= curratios.length;



		    if ((curratios.length > 0)&&(prevratios.length > 0))
		    {
	               dpval = theMannWhitneyUTest.mannWhitneyUTest(curratios,prevratios);
		    }
		    else
		    {
		       dpval = 1;
		    }

		    Arrays.sort(curratios);

		    if (curratios.length == 0)
		    {
			medians[nfile][nseq][npos] = 0;
		    }
		    else if (curratios.length % 2 == 1)
		    {
			medians[nfile][nseq][npos] = curratios[curratios.length/2];
		    }
		    else
		    {
			medians[nfile][nseq][npos] = (curratios[curratios.length/2]+curratios[curratios.length/2-1])/2;
		    }

		    //System.out.println("-->"+nfile+"\t"+nseq+"\t"+npos+"\t"+dpval+"\t"+davgprev+"\t"+davgcurr);//+pvals[nfile].length+"\t");

		    if ((prevratios.length == 0)||(curratios.length == 0))
		    {
			pvalsGREATER[nfile][nseq][npos]= 1;
			pvalsLESS[nfile][nseq][npos]= 1;
		    }
		    else if (davgprev >= davgcurr)
		    {
		       pvalsLESS[nfile][nseq][npos] = dpval/2;
		       pvalsGREATER[nfile][nseq][npos] = 1-dpval/2.0;
		    }
		    else
		    {
		       pvalsGREATER[nfile][nseq][npos] = dpval/2;
		       pvalsLESS[nfile][nseq][npos] = 1-dpval/2.0;
		    }
		    //System.out.println(nfile+"\t"+nseq+"\t"+npos+"\t"+dpval+"\t"+curratios.length+"\t"+prevratios.length);		    
		    prevratios = curratios;
		}
	    }
	}      
    
	//combinedpvals = new double[pvals[0].length][pvals[0][0].length];

	for (int na = 0; na < pvalsGREATER[0].length; na++)
	{
	    for (int nb = 1; nb < pvalsGREATER[0][na].length; nb++)
	    {
		//System.out.print(na+"\t"+nb);
		double dsum = 0;

		for (int nc = 0; nc < pvalsGREATER.length; nc++)
		{
		   dsum += Math.log(pvalsGREATER[nc][na][nb]);
		}

		double dtest = -2*dsum;
		ChiSquaredDistribution theChiSq = new ChiSquaredDistribution(2*pvalsGREATER.length);
		double dpval = 1- theChiSq.cumulativeProbability(dtest);
		combinedpvalsGREATER[na][nb] = dpval;

                dsum = 0;
                for (int nc = 0; nc < pvalsLESS.length; nc++)
	        {
		    // System.out.println(nc+"\t"+na+"\t"+nb+"\t"+pvalsLESS[nc][na][nb]);
	      	   dsum += Math.log(pvalsLESS[nc][na][nb]);
		}
                dtest = -2*dsum;
                theChiSq = new ChiSquaredDistribution(2*pvalsLESS.length);
                dpval = 1- theChiSq.cumulativeProbability(dtest);
                combinedpvalsLESS[na][nb] = dpval;
		//multiply by two is for multiple testing correction

	    }
	 }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////


    static class RankRec
    {
       double dval;
       boolean bprev;

       RankRec(double dval, boolean bprev)
       {
	  this.dval = dval;
 	  this.bprev = bprev;
       }
    }

    static class RankRecCompare implements Comparator
    {

       public int compare(Object o1, Object o2)
       {
	  RankRec r1 = (RankRec) o1;
	  RankRec r2 = (RankRec) o2;

	  if (r1.dval < r2.dval)
	  {
	     return -1;
 	  }
	  else if (r1.dval > r2.dval)
	  {
	     return 1;
	  }
	  else
	  {
	     return 0;
	  }
       }
    }


    /*
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    private void createrandomized(ArrayList[][][] alvals, ArrayList[][][] randomvals, Random theRandom)
    {
	ArrayList aldoublevals = new ArrayList();
	for (int na = 0; na < alvals.length; na++)
	{
	    for (int nb = 0; nb < alvals[na].length; nb++)
	    {
		for (int nc = 0; nc < alvals[na][nb].length; nc++)
	        {
		    randomvals[na][nb][nc] = new ArrayList();
		    for (int nd = 0; nd < alvals[na][nb][nc].size(); nd++)
		    {
			aldoublevals.add(new RandomRec(((Double) alvals[na][nb][nc].get(nd)).doubleValue(),theRandom.nextDouble()));
		    }
		}
	    }
	}

	RandomRec[] aldoublevalsA = new RandomRec[aldoublevals.size()];
	for (int nk = 0; nk < aldoublevalsA.length; nk++)
	{
	    aldoublevalsA[nk] = (RandomRec) aldoublevals.get(nk);
	}

	Arrays.sort(aldoublevalsA, new RandomRecCompare());
	int ne = 0;
	for (int na = 0; na < alvals.length; na++)
	{
	    for (int nb = 0; nb < alvals[na].length; nb++)
	    {
		for (int nc = 0; nc < alvals[na][nb].length; nc++)
		{
		    for (int nd = 0; nd < alvals[na][nb][nc].size(); nd++)
		    {
			randomvals[na][nb][nc].add(((RandomRec) (aldoublevalsA[ne])).dval);
			//System.out.println(((RandomRec) (aldoublevalsA[ne])).dval);
			ne++;
		    }
		}
	    }
	}
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static class RandomRec
    {
	double drandomval;
	double dval;

	public RandomRec(double dval,double drandomval)
	{
	    this.drandomval = drandomval;
	    this.dval = dval;
	}
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    static class RandomRecCompare implements Comparator
    {

	public int compare(Object o1, Object o2)
	{

	    RandomRec r1 = (RandomRec) o1;
	    RandomRec r2 = (RandomRec) o2;

	    if (r1.drandomval < r2.drandomval)
	    {
		return -1;
	    }
	    else if (r1.drandomval > r2.drandomval)
	    {
		return 1;
	    }
	    else
	    {
		return 0;
	    }
	}
    }

    */

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    public static void main(String[] args) throws IOException
    {


        int nargindex = 0;


	boolean bok = true;

	String szcommand = "";
	if (args.length >=1)
	{
           szcommand = args[nargindex++];
	}

        if (szcommand.equalsIgnoreCase("Version"))
	{
	   System.out.println("This is version 1.0.1 of SHARPR");
	}
	else if (szcommand.equals("ExecuteAll"))
	{
	    String szrnafilelist = null;
            String szdnafilelist = null;
            String szinterpolateoutputfile = null;

	    int ndnacutoff = SHARPR.DEFAULT_THRESHOLD;
	    int nrnapseudocount = SHARPR.DEFAULT_RNAPSEUDO;
	    int ndnapseudocount = SHARPR.DEFAULT_DNAPSEUDO;
	    int ntileindex = 0;
	    double dvarprior1 = SHARPR.DEFAULT_VARPRIOR1;// 0;
            double dvarprior2 = SHARPR.DEFAULT_VARPRIOR2;// 0;
	    int ntilelength = 0;
	    int nstepsize = 0;
	    int numtiles = 0;


	   while (nargindex < args.length-7)
	   {
              //imputation file is CELL type, then mark, then data file
	      //default output is the full imputation grid
	      //the -r option specifies the output resolution to be used

	      //CELL type
	      //files should be organized CELL type then directory

	      //input is a cell type/mark text file indicating
	      //if that cell type and mark should also
	      String szoption;

	      szoption = args[nargindex++];


	      if (szoption.equals("-c"))
	      {
	         ndnacutoff = Integer.parseInt(args[nargindex++]);
	      }
	      else if (szoption.equals("-p"))
	      {
		 nrnapseudocount = Integer.parseInt(args[nargindex++]);
	      }
	      else if (szoption.equals("-d"))
	      {
		  ndnapseudocount = Integer.parseInt(args[nargindex++]);
	      }
              else if (szoption.equals("-v1"))
	      {
		  dvarprior1 = Double.parseDouble(args[nargindex++]);
	      }
              else if (szoption.equals("-v2"))
	      {
		  dvarprior2 = Double.parseDouble(args[nargindex++]);
	      }
	      else
	      {
	         bok = false;
	      }
	   }

           if (nargindex == args.length-7)
	   {
              //input data
	      szrnafilelist = args[nargindex++];
	      szdnafilelist = args[nargindex++];
	      ntileindex = Integer.parseInt(args[nargindex++]);
	      //dvarprior1 = Double.parseDouble(args[nargindex++]);
              //dvarprior2 = Double.parseDouble(args[nargindex++]);
	      ntilelength = Integer.parseInt(args[nargindex++]);
	      nstepsize = Integer.parseInt(args[nargindex++]);
	      numtiles = Integer.parseInt(args[nargindex++]);
	      szinterpolateoutputfile = args[nargindex++];
	   }
           else
           {
	       bok = false;
	   }

           if (!bok)
	   {
	       System.out.println("USAGE: java ExecuteAll [-c cutoff][-p pseudcountrna][-d pseudocountdna][-v1 varprior1][-v2 varprior2] rnafilelist dnafilelist npositionindex "+
                                   "tilelength stepsize numtilepos interpolateoutfile");

	   }
	   else
	   {

              try
	      {
	          new SHARPR(szrnafilelist, szdnafilelist, ntileindex, dvarprior1, dvarprior2, ntilelength, nstepsize, numtiles, 
			     szinterpolateoutputfile,ndnacutoff,nrnapseudocount,ndnapseudocount);
	      }
              catch (IOException ioex)
              {
      	         ioex.printStackTrace(System.out);
	      }
	   }

	}
        else if (szcommand.equalsIgnoreCase("Normalize"))
	{

	   String szrnafile = null;
	   String szdnafile = null;

	   int ndnacutoff = SHARPR.DEFAULT_THRESHOLD;
	   int nrnapseudocount = SHARPR.DEFAULT_RNAPSEUDO;
	   int ndnapseudocount = SHARPR.DEFAULT_DNAPSEUDO;
	   //int nnormalizewidth = SHARPR.DEFAULT_NORMALIZEWIDTH;

	   String sznormoutfile = null;
	   String sznormoutfilenoavg = null;



	   while (nargindex < args.length-3)
	   {
              //imputation file is CELL type, then mark, then data file
	      //default output is the full imputation grid
	      //the -r option specifies the output resolution to be used

	      //CELL type
	      //files should be organized CELL type then directory

	      //input is a cell type/mark text file indicating
	      //if that cell type and mark should also
	      String szoption;

	      szoption = args[nargindex++];

	      if (szoption.equals("-c"))
	      {
	         ndnacutoff = Integer.parseInt(args[nargindex++]);
	      }
	      else if (szoption.equals("-p"))
	      {
		 nrnapseudocount = Integer.parseInt(args[nargindex++]);
	      }
	      else if (szoption.equals("-n"))
	      {
		  sznormoutfilenoavg = args[nargindex++];
	      }
	      else if (szoption.equals("-d"))
	      {
		  ndnapseudocount = Integer.parseInt(args[nargindex++]);
	      }
	      else
	      {
	         bok = false;
	      }
	   }

           if (nargindex == args.length-3)
	   {
              //input data
	      szrnafile = args[nargindex++];
	      szdnafile = args[nargindex++];
	      sznormoutfile = args[nargindex++];
	   }
           else
           {
	       bok = false;
	   }

	   //System.out.println(sznormoutfilenoavg);
           if (!bok)
	   {
     	      System.out.println("USAGE: java SHARPR Normalize [-c cutoff][-p pseudcountrna][-d pseudocountdna][-n noavgoutfile] rnafile dnafile outfile");
	   }
	   else
	   {

              try
	      {
	          new SHARPR(szrnafile, szdnafile, sznormoutfile, sznormoutfilenoavg, ndnacutoff, nrnapseudocount,ndnapseudocount);
	      }
              catch (IOException ioex)
              {
      	         ioex.printStackTrace(System.out);
	      }
	   }
	} 
        else if ((szcommand.equalsIgnoreCase("Deconvolve"))||(szcommand.equalsIgnoreCase("Infer")))
	{
	    boolean bstandardize = true;
	    String szdeconvolveinputfile = null;
	    String szdeconvolveoutputfile = null;
	    double dpriorvar = 0;
	    int ntilelength = 0;
	    int nstepsize = 0;
	    int numtiles = 0;

	   while (nargindex < args.length-6)
	   {

	      String szoption;

	      szoption = args[nargindex++];

	      if (szoption.equals("-nostandardize"))
	      {
		  bstandardize = false;
	      }
	      else
	      {
	         bok = false;
	      }
	   }

	   if (nargindex == args.length - 6)
	   {
	      szdeconvolveinputfile = args[nargindex++];
	      szdeconvolveoutputfile = args[nargindex++];
	      dpriorvar = Double.parseDouble(args[nargindex++]);
	      ntilelength = Integer.parseInt(args[nargindex++]);
	      nstepsize = Integer.parseInt(args[nargindex++]);
	      numtiles = Integer.parseInt(args[nargindex++]);
	   }
	   else
	   {
	      bok = false;
	   }

           if (!bok)
	   {
     	      System.out.println("USAGE: java SHARPR Infer [-nostandardize] inputtablefile outputfile varprior tilelength stepsize numtilepos");
	   }
	   else
	   {
              try
	      {
		  //    new SHARPR(szdeconvolveinputfile, szdeconvolveoutputfile, dpriorvar, ntilelength, nstepsize, numtiles, ngibbsburnin, ngibbstotaliterations, bstandardize);
	          new SHARPR(szdeconvolveinputfile, szdeconvolveoutputfile, dpriorvar, ntilelength, nstepsize, numtiles, bstandardize);
	      }
              catch (IOException ioex)
              {
      	         ioex.printStackTrace(System.out);
	      }
	   }
	}
	else if (szcommand.equals("Enrich"))
	{

	    double dvaluebins = SHARPR.DEFAULT_VALUEBINS;
	    int nbasetooverlap = -1;
	    int nmincountextreme =  SHARPR.DEFAULT_MINCOUNTEXTREME;
	    int numoverlapbins = SHARPR.DEFAULT_NUMOVERLAPBINS;
	    String szreportercoordinates = null;
	    String szoverlapcoordinates = null;
	    String szbasepredictions = null;
	    String  szoverlapoutput = null;
	    String szfeatureoutput = null;
	    String szsummaryoutput = null;
	    boolean bmaxenrich = false;
	    boolean bvaluebins = false;
	    boolean bneeddefaultbins = true;

	   while (nargindex < args.length-4)
	   {
	       String szoption;

	       szoption = args[nargindex++];

	       if (szoption.equals("-f"))
	       {
		   szfeatureoutput = args[nargindex++];
	       }
	       else if (szoption.equals("-b"))
	       {
		   nbasetooverlap = Integer.parseInt(args[nargindex++]);
	       }
	       else if (szoption.equals("-maxabsenrich"))
	       {
		   bmaxenrich = true;
	       }
	       else if (szoption.equals("-s"))
	       {
		   szsummaryoutput = args[nargindex++];
	       }
	       else if (szoption.equals("-v"))
	       {
		   dvaluebins = Double.parseDouble(args[nargindex++]);
	       }
	       else if (szoption.equals("-n"))
	       {
		   numoverlapbins = Integer.parseInt(args[nargindex++]);
		   bneeddefaultbins = false;
	       }
	       else if (szoption.equals("-valuebins"))
	       {
		   bvaluebins = true;
	       }
	       else if (szoption.equals("-e"))
	       {
                   nmincountextreme = Integer.parseInt(args[nargindex++]);
	       }
	       else
	       {
		   bok = false;
	       }
	   }

	   if ((bneeddefaultbins)&& (bmaxenrich))
	   {
	       numoverlapbins = SHARPR.DEFAULT_NUMOVERLAPBINS_MAX;
	   }

	    if (nargindex == args.length -4)
	    {
		szreportercoordinates = args[nargindex++];
		szoverlapcoordinates = args[nargindex++];
		szbasepredictions = args[nargindex++];
		szoverlapoutput = args[nargindex++];
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("USAGE: java SHARPR Enrich [-b basetooverlap][-e mincountextreme][-f featureoutput][-maxabsenrich][-n numoverlapbins][-s summaryoutput][-valuebins][-v value] reportercoordinates overlapcoordinates basepredictions overlapoutput"); 
                //Combine fileset1 [fileset2] outputfile");
	    }
	    else
	    {
                try
		{
		    new SHARPR(nbasetooverlap, numoverlapbins, szreportercoordinates, szoverlapcoordinates, szbasepredictions, szoverlapoutput, szfeatureoutput, 
                                  bmaxenrich,szsummaryoutput, dvaluebins, bvaluebins,nmincountextreme);
		}
		catch (IOException ioex)
	        {
	      	   ioex.printStackTrace(System.out);
	        }
	    }
	}
	else if (szcommand.equalsIgnoreCase("Combine"))
	{

	    String szcombinefileset1 = null;
	    String szcombinefileset2 = null;
	    String szoutputcombine = null;

	    while (nargindex < args.length-2)
	    {
		String szoption;

		szoption = args[nargindex++];

		if (szoption.equals("-c"))
		{
		    szcombinefileset2 = args[nargindex++];
		}
		else
		{
		    bok = false;
		}
	    }

	   if (nargindex == args.length -2)
	   {
               szcombinefileset1 = args[nargindex++];
	       szoutputcombine = args[nargindex++];
	   }
	   /*
           //changed here so second file set is an option not specified at the command line
	   else if (nargindex == args.length - 3)
	   {
	       szcombinefileset1 = args[nargindex++];
	       szcombinefileset2 = args[nargindex++];
	       szoutputcombine = args[nargindex++];
	   }
	   */
	   else
	   {
	       bok = false;
	   }


           if (!bok)
	   {
	      System.out.println("USAGE: java Combine [-c fileset2] fileset1 outputfile");
	   }
           else
	   {
              try
	      {
		  new SHARPR(szcombinefileset1, szcombinefileset2, szoutputcombine);
	      }
              catch (IOException ioex)
	      {
		  ioex.printStackTrace(System.out);
	      }
	   }
	}
	else if (szcommand.equalsIgnoreCase("Interpolate"))
	{

	    String szinterpolateinputfile = null;
	    String szinterpolateoutputfile = null;
	    int nstepsize = 0;

	   if (nargindex == args.length-3)
	   {
              //input data
	      //String szconvertinputfile, String szconvertoutputfile
	       szinterpolateinputfile = args[nargindex++];
	       szinterpolateoutputfile = args[nargindex++];
	       nstepsize = Integer.parseInt(args[nargindex++]);
	   }
	   else
	   {
              bok = false;
	   }


           if (!bok)
	   {
	      System.out.println("USAGE: java SHARPR Interpolate interpolateinputfile interpolateoutputfile stepsize");
	   }
           else
	   {
              try
	      {
		  new SHARPR(nstepsize,szinterpolateinputfile, szinterpolateoutputfile);
	      }
              catch (IOException ioex)
	      {
		  ioex.printStackTrace(System.out);
	      }
	   }
	}
	else if (szcommand.equalsIgnoreCase("AdjacentChanges"))
	{
	    String szdatafiles = null;
	    int numtiles = 0;
	    int nextension = SHARPR.DEFAULT_TILEEXTENSION;
	    int nstepsize = 0;
	    String szsequences = null;
	    String szreportercoordinates = null;
	    String szadjacentoutfile = null;

	    while (nargindex < args.length-6)
	    {
	       String szoption;

	       szoption = args[nargindex++];

	       if (szoption.equals("-e"))
	       {
	          nextension = Integer.parseInt(args[nargindex++]);
	       }
	       else
	       {
	          bok = false;
	       }
	    }


	    if (nargindex == args.length - 6)
	    {
		szdatafiles = args[nargindex++];
		nstepsize = Integer.parseInt(args[nargindex++]);
		numtiles = Integer.parseInt(args[nargindex++]);
		szsequences = args[nargindex++];
		szreportercoordinates= args[nargindex++];
		szadjacentoutfile = args[nargindex++];
	    }
	    else
	    {
		bok = false;
	    }

	    if (!bok)
	    {
		System.out.println("USAGE: java AdjacentChanges [-e extension] datafiles stepsize numpositions sequences coordinates outfile");
	    }
	    else
	    {
               try
	       {
		   new SHARPR(szdatafiles, nstepsize, numtiles, szsequences, nextension, szreportercoordinates, szadjacentoutfile);
	       }
	       catch (IOException ioex)
	       {
	          ioex.printStackTrace(System.out);
	       }
	    }
	}
	else if (szcommand.equalsIgnoreCase("ConvertTable"))
	{
	    int ntileindex = 0;
	    String szconvertinputfile = null;
	    String szconvertoutputfile = null;

	    boolean brecentermean = false;
	    boolean brecenterends = false;


	   while (nargindex < args.length-3)
	   {

	      String szoption;

	      szoption = args[nargindex++];

	      if (szoption.equals("-recentermean"))
	      {
		  brecentermean = true;
	      }
	      else if (szoption.equals("-recenterends"))
	      {
		  brecenterends = true;
	      }
	      else
      	      {
       	          bok = false;
	      }

	   }

	   if (nargindex == args.length-3)
	   {
	      //input data
	       //String szconvertinputfile, String szconvertoutputfile
	       szconvertinputfile = args[nargindex++];
	       szconvertoutputfile = args[nargindex++];
	       ntileindex = Integer.parseInt(args[nargindex++]);
	   }
	   else
	   {
	       bok = false;
	   }


	   if (!bok)
	   {
	       System.out.println("USAGE: java SHARPR ConvertTable [-recentermean][-recenterends] convertinputfile convertoutputfile npositionindex");
	   }
	   else
	   {
              try
	      {
		  new SHARPR(szconvertinputfile, szconvertoutputfile, ntileindex, brecentermean, brecenterends);
	      }
              catch (IOException ioex)
	      {
	         ioex.printStackTrace(System.out);
	      }
	   }
	}
	else 
	{
	    System.out.println("Need to specify mode ExecuteAll|Normalize|ConvertTable|Infer|Combine|Interpolate|Enrich|AdjacentChanges|Version");
	}

    }

}
