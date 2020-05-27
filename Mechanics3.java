package research;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class Mechanics3 {
		
		//constants
		private static int HR = 60; //bpm
		private static double tc = 60.0/HR; //amt sec in 1 cardiac cycle
		private static double Tmax = 0.2 + 0.15*tc; //time when Emax occurs
		private static double Rs = 1.0;    //mmHg/mL systemic vascular resistance
		private static double Rm = 0.0050; //mmHg/mL mitral valve resistance
		private static double Ra = 0.0010; //mmHg/mL aortic valve resistance
		private static double Rc = 0.0398; //mmHg/mL characteristic resistance
		private static double Cr = 4.4000; //mL/mmHg left atrial compliance
		private static double Cs = 1.3300; //mL/mmHg systemic compliance
		private static double Ca = 0.0800; //mL/mmHg aortic compliance
		private static double Ls = 0.0005; //mmHg^2/mL inertance of blood in aorta */
		private static double Emax = 2;    //max elastance
		private static double Emin = 0.06; //min elastance
		private static double V0 = 10.0;   //reference vol of blood under 0 pressure 
		
	public static void main(String[] args) throws IOException {
		
		PrintWriter pw = new PrintWriter(new FileWriter("all_data.out"));
		PrintWriter pw1 = new PrintWriter(new FileWriter("LVP.out"));
		PrintWriter pw2 = new PrintWriter(new FileWriter("LAP.out"));
		PrintWriter pw3 = new PrintWriter(new FileWriter("AP.out"));
		PrintWriter pw4 = new PrintWriter(new FileWriter("AOP.out"));
		PrintWriter pw5 = new PrintWriter(new FileWriter("Qa.out"));
		PrintWriter pw6 = new PrintWriter(new FileWriter("LVV.out")); 
		
		//entering user-input values
		Scanner keyboard = new Scanner(System.in);
		System.out.println("Enter the heart rate (in bpm). Enter -1 if using default");
		HR = keyboard.nextInt();
		if(HR==-1)
			HR = 60;
		System.out.println("Enter the number of cycles");
		tc = (double) keyboard.nextInt()*60/HR;
		System.out.println("Enter the Systemic Vascular Resistance (in mmHg*s/mL). Enter -1 if using default");
		Rs = keyboard.nextDouble();
		if(Rs==-1)
			Rs = 1.0;
		System.out.println("Enter the Mitral Valve Resistance (in mmHg*s/mL). Enter -1 if using default");
		Rm = keyboard.nextDouble();
		if(Rm==-1)
			Rm = 0.005;
		System.out.println("Enter the Aortic Valve Resistance (in mmHg*s/mL). Enter -1 if using default");
		Ra = keyboard.nextDouble();
		if(Ra==-1)
			Ra = 0.001;
		System.out.println("Enter the Characteristic Resistance (in mmHg*s/mL). Enter -1 if using default");
		Rc = keyboard.nextDouble();
		if(Rc==-1)
			Rc = 0.0398;
		System.out.println("Enter the Left Atrial Compliance (in mL/(mmHg*s)). Enter -1 if using default");
		double Cr = keyboard.nextDouble();
		if(Cr==-1)
			Cr = 4.4000;
		System.out.println("Enter the Systemic Compliance (in mL/(mmHg*s)). Enter -1 if using default");
		Cs = keyboard.nextDouble();
		if(Cs==-1)
			Cs = 1.3300;
		System.out.println("Enter the Aortic Compliance (in mL/(mmHg*s)). Enter -1 if using default");
		Ca = keyboard.nextDouble();
		if(Ca==-1)
			Ca = 0.0800;
		System.out.println("Enter the Aortic Blood Inertance (in mmHg^2/mL). Enter -1 if using default");
		Ls = keyboard.nextDouble();
		if(Ls==-1)
			Ls = 0.0005; 
		
		double Tmax = 0.2 + 0.15*tc;		
		
		double dt = 0.0001; //discrete time step
		
		//initial values
		double[] LVP = new double[10002]; //left ventricular pressure
		double[] LAP = new double[10002]; //left atrial pressure
		double[] AP = new double[10002]; //arterial pressure
		double[] AOP = new double[10002]; //aortic pressure
		double[] Qa = new double[10002]; //flow rate
		double[] Cv = new double[10002]; //compliance
		
		System.out.println("Enter the initial Left Ventricular Pressure (in mmHg). Enter -1 if using default");
		LVP[0] = keyboard.nextDouble();
		if(LVP[0]==-1)
			LVP[0] = 8.2;
		
		System.out.println("Enter initial Left Atrial Pressure (in mmHg). Enter -1 if using default");
		LAP[0] = keyboard.nextDouble();
		if(LAP[0]==-1)
			LAP[0] = 7.6;
		
		System.out.println("Enter initial Arterial Pressure (in mmHg). Enter -1 if using default");
		AP[0] = keyboard.nextDouble();
		if(AP[0]==-1)
			AP[0] = 67.0;
		
		System.out.println("Enter initial Aortic Pressure (in mmHg). Enter -1 if using default");
		AOP[0] = keyboard.nextDouble();
		if(AOP[0]==-1)
			AOP[0] = 80; //default value
		
		System.out.println("Enter initial Flow Rate (in mL/s). Enter -1 if using default");
		Qa[0] = keyboard.nextDouble();
		if(Qa[0]==-1)
			Qa[0] = 8.2; //default value
		
		System.out.println("Enter initial Ventricular Compliance (in mL/mmHg). Enter -1 if using default");
		Cv[0] = keyboard.nextDouble();
		if(Cv[0]==-1)
			Cv[0] = 8.2; //default value
		
		double t = 0;
		double[][] dx = new double[5][1];
		
		int i = -1;
		double tn = 0.0; //proportion of total time that has passed (t/Tmax)
		double En = 0.0; //normalized elastance
		double Dm = 0.0; //mitral diode on or off
		double Da = 0.0; //aortic dioide on or off
		double dCv = 0; //differentiated compliance
		double[] Tn = new double[10002]; //array of tn
		Tn[0] = t;
		double[] E = new double[10002]; //array for elastance over time
		double[] LVV = new double[10002]; //left ventricular volume over time
		double[] Time = new double[10002]; //counting the discrete time steps; for graphing purposes
		
		//matrices for eventual eventual derivative calculation, values not final and will be altered in while loop
		double[][] A = {
				{0,0,0,0,0},
		        {0,-1/(Rs*Cr),1/(Rs*Cr),0,0},
		        {0,1/(Rs*Cs),-1/(Rs*Cs),0,1/Cs},
		        {0,0,0,0,-1/Ca},
		        {0,0,-1/Ls,1/Ls,-Rc/Ls}
		};
		double[][] B = {
				{0,0},
				{-1/Cr,0},
				{0,0},
				{0,1/Ca},
				{0,0}
		};
		double[][] C = new double[2][1];
		
		while(t<=tc) {
			i++;
			tn = (t-Math.floor(t))/Tmax;
			
			//normalized elastance from double hill equation
			En = 1.55*Math.pow(tn/0.7, 1.9)/(1+Math.pow(tn/0.7, 1.9))*(1/(1+Math.pow(tn/1.17, 21.9)));
			
			E[i] = (Emax-Emin)*En + Emin;
			LVV[i] = LVP[i]/E[i] + V0;
			Cv[i] = 1/E[i];
			
			if(LAP[i] > LVP[i]) //determine if mitral valve opens
				Dm = 1;
			else
				Dm = 0;
			if(LVP[i] > AOP[i]) //determine if aortic valve opens
				Da = 1;
			else
				Da = 0;
			if(i>1) 
				dCv = (Cv[i]-Cv[i-1])/dt; //change in compliance
			else
				dCv = 0;
			
			A[0][0] = -dCv/Cv[i];
			B[0][0] = 1/Cv[i];
			B[0][1] = -1/Cv[i];
			C[0][0] = (Dm/Rm)*(LAP[i] - LVP[i]);
			C[1][0] = (Da/Ra)*(LVP[i] - AOP[i]);

			
			copyMatrix(dx, addMatrix(multiplyMatrix(A, makeMatrix(LVP,LAP,AP,AOP,Qa,i)),multiplyMatrix(B,C)));
			multiplyMatrixByConstant(dx,dt);
			LVP[i+1] = LVP[i] + dx[0][0];
			LAP[i+1] = LAP[i] + dx[1][0];
			AP[i+1] = AP[i] + dx[2][0];
			AOP[i+1] = AOP[i] + dx[3][0];
			Qa[i+1] = Qa[i] + dx[4][0];
			
			t += dt;
			Time[i+1] = t;
		}
		
		i++;
		tn = (t-Math.floor(t))/Tmax;
		En = 1.55*Math.pow(tn/0.7, 1.9)/(1+Math.pow(tn/0.7, 1.9))*(1/(1+Math.pow(tn/1.17, 21.9)));
		E[i] = (Emax-Emin)*En + Emin;
		LVV[i] = LVP[i]/E[i] + V0;
		
		for(int n=0; n<Time.length; n++) {
			pw.println(Time[n]+"\t"+LVP[n]+"\t"+LAP[n]+"\t"+AP[n]+"\t"+AOP[n]+"\t"+Qa[n]+"\t"+LVV[n]);
		} 
		
		for(int n=0; n<Time.length; n++)
			pw1.println(Time[n]+"\t"+LVP[n]);
		
		for(int n=0; n<Time.length; n++) {
			pw2.println(Time[n]+"\t"+LAP[n]);
		}
		for(int n=0; n<Time.length; n++) {
			pw3.println(Time[n]+"\t"+AP[n]);
		}
		for(int n=0; n<Time.length; n++) {
			pw4.println(Time[n]+"\t"+AOP[n]);
		}
		for(int n=0; n<Time.length; n++) {
			pw5.println(Time[n]+"\t"+Qa[n]);
		}
		for(int n=0; n<LVP.length; n++) {
			pw1.println(LVP[n]+"\t"+LVV[n]);
		} 
		
		pw.close();
		
		pw2.close();
		pw3.close();
		pw4.close();
		pw5.close();
		pw6.close(); 
	}
	
	public static double[][] multiplyMatrix(double[][] firstMatrix, double[][] secondMatrix) {
        int r1 = firstMatrix.length;
        int c1 = firstMatrix[0].length;
        int r2 = secondMatrix.length;
        int c2 = secondMatrix[0].length;
		double[][] product = new double[r1][c2];
        for(int i = 0; i < r1; i++) {
            for (int j = 0; j < c2; j++) {
                for (int k = 0; k < c1; k++) {
                    product[i][j] += firstMatrix[i][k] * secondMatrix[k][j];
                }
            }
        }

        return product;
    }
	
	public static void multiplyMatrixByConstant(double[][] a, double con) {
		for(int r=0; r<a.length; r++)
			for(int c=0; c<a[0].length; c++)
				a[r][c] = a[r][c]*con;
	}
	
	public static double[][] addMatrix(double[][] a, double[][] b) {
		double[][] result = new double[a.length][a[0].length];
		for(int r=0; r<a.length; r++)
			for(int c=0; c<a[0].length; c++)
				result[r][c] = a[r][c]+b[r][c];
		return result;
	}
	
	public static double[][] makeMatrix(double[] a, double[] b, double[] c, double[] d, double[] e, int index) {
		double[][] result = new double[5][1];
		result[0][0] = a[index];
		result[1][0] = b[index];
		result[2][0] = c[index];
		result[3][0] = d[index];
		result[4][0] = e[index];
		return result;
	}
	
	public static void copyMatrix(double[][] original, double[][] template) {
		for(int r=0; r<original.length; r++)
			for(int c=0; c<original[0].length; c++)
				original[r][c] = template[r][c];
	}
}
