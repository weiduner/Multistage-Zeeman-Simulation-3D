
//Zeeman deceleration for oxygen molecules by using Sida's output 3d magnetic field
//Molecules climb the potential hill from zero phase angle instead of -90 degree to compare with Sida's simulation at the same condition

//simulation for experimental data of April 29th

#include "stdlib.h"
#include "math.h"
#include "time.h"
#include "string.h"
#include "stdio.h"
#include "fstream.h"
#include "iomanip.h"
#include "iostream.h"

void MoleculeEnsemble(double *MoleculeList);

double PartialEnergyToB(double n);

void MToAcceleration(double coeff,double *GradB,double x,double y,double z,double N,double *ax, double *ay, double *az);
int Interpolation(double *a,double dx,double dy,double dz, double *b);
void Runge_Kutta(double coeff,double *GradB,double x0,double y0,double z0,double vx0,double vy0,double vz0,double N,double *x,double *y,double *z,double *vx,double *vy,double *vz,double dt);

double DerivativeInterpolation(double *DiffB,double f,double N,double J,double Mj);

void FreeFly(double Time[],int ToFSignal[],char *ptwriteZeemanToF);
void FiniteSwitchingSequence(double *Bn,double *DiffB,double *GradB,double TimeArray[],double VExit[],double *ta,char* ptwriteTimeSequence);

double SwitchingField(double t,double t0);
double ReadTiming(char* PRead,double a[]);
void end_on(double t,double z);
void end_off(double t);
void save_initial();
void trace(double TimeArray[],double *Bn);
void save_end(double t);

void NewFiniteSwitchingSequence(double *GradB,double *fitting,double *Bx_PM,double *By_PM,double *Bz_PM,double TimeOn[],double TimeArray[],double VExit[],double *ta,char* ptwriteTimeSequence);


void save_TOF(double t);
void New_SlowerFlight(double TimeOn[],double TimeArray[],double Time[],int ToFSignal[],double *GradB,double *fitting,char *ptwriteZeemanToF);
double NewFittingSwitchingField(double *fitting,int i,double t,double ton);
double AppFittingSwitchingField(double *fitting,int i,double t,double ton,double toff);
double ReadTiming(char* PRead,double a[],double b[],double c[]);

void Trap_Runge_Kutta(double *Btx,double *Bty,double *Btz,double x0,double y0,double z0,double vx0,double vy0,double vz0,double m,double *x,double *y,double *z,double *vx,double *vy,double *vz,double dt);
int InterpolationX_trap(double *a,double dx,double dy,double dz, double *b,double *c);
int InterpolationY_trap(double *a,double dx,double dy,double dz, double *b,double *c);
int InterpolationZ_trap(double *a,double dx,double dy,double dz, double *b,double *c);
void TrapAcceleration(double *Btx,double *Bty,double *Btz,double x,double y,double z,double m,double *ax, double *ay, double *az);


const double KB=1.380650524e-23;//Boltzmann constant 
const double h         =       6.626069311e-34;
const double c			=		299792458;
const double PI				=	3.14159265358979323846264338328;
const double amu=1.66053886e-27;//atomic mass unit
//const double DebyeToColoum = 3.33564e-30;


clock_t start_time,elapsed;

double elapsed_time;

double *Bz=NULL; //magnetic field along z direction
double *Br=NULL; //magnetic field along r direction
double *Bn=NULL; //total magnetic field
double *Shift=NULL;  //gradient of magnetic field
double *DiffB=NULL;  //magnetic field along z direction
double *GradB=NULL;  //gradient of magnetic field
//double *GradBx=NULL;  //gradient of magnetic field along x direction
//double *GradBy=NULL;  //gradient of magnetic field along y direction
//double *GradBz=NULL;  //gradient of magnetic field along z direction

//Field data from Mario
double ZMin=0;
double ZMax=12e-3;
double RMin=-2e-3;
double RMax=2e-3;
int ZNumber=501;//int((ZMax-ZMin)/ZGridSize+1)
int RNumber=51;//int((RMax-RMin)/RGridSize+1)
double RGridSize = 8e-5;
double ZGridSize = 2.4e-5;


//for the permanent trap
double GridSize0 = 0.1e-3;
double XMin0=-5e-3;
double XMax0=5e-3;
double YMin0=-5e-3;
double YMax0=5e-3;
double ZMin0=-15e-3;
double ZMax0=15e-3;
int XNumber0=101;
int YNumber0=101;
int ZNumber0=301;	//3070501


//double field = 4.27;			//for 600A and ~68us current pulse
double field = 4.27;			//for ~80us current pulse 
//FieldScaleFactor is the factor for scaling the maximum field due to inperfection of calculation for the magnetic spatial distribution
double FieldScaleFactor = 1.037*105/140;
//double FieldScaleFactor = 1.037*180/156;

double FieldStep = 1e-4; //1 Gauss
double FieldMax = 6;	//6 Tesla
double MinScale = 1;		//corresponding to 550A
double MaxScale = 1.3;		//corresponding to 750A
int NumOfField = (int)(floor)(FieldMax/FieldStep);


//in experiments we always meet the situation that the maximum field strength is different from what we have calculated from
//Finite Element Software, then we have to scale it to the one we have calculted.
//Esp. for Zeeman deceleration, differences always exist between different measurements or between measurements and calculations

//for calculation of Zeeman shift of the specific quantum state
double mSample = 26; 
double mass = mSample*amu;	

double B = 43100.46e6;//rotational constant for oxygen in the 3Sigma_g, v=0 ground state
//////////////////////////

/******Distribution Parameters *******/


long MoleculeNumber = long(5e6 );

double Xcenter = 0;
double Ycenter = 0;
double Zcenter = 2e-3;	//distance from discharge plate to exit of the nozzle

double DeltaX = 33e-6;     //the nozzle aperture is 200um 
double DeltaY = 33e-6;
double DeltaZ = 0.005;      //(opening time of nozzle * central speed)/6


int FlatX = 0;
int FlatY = 0;
int FlatZ = 0;

//0.1 is about the same with 0.3m/s in terms of signal intensity of 43m/s
double DeltavX = 0.1; //calculated from the geometry: distance from nozzle to skimmer (to 1st coil), diameter of skimmer, inner radius of Zeeman coil
double DeltavY = 0.1;        
double DeltavZ = 38; 

double vXcenter = 0;
double vYcenter = 0;
double vZcenter = 400;

double Vsyn = 960;

double Phi = 10.1;	//used for first N stages

/*******Overtone Parameters********/

int NumberOfStages = 85; 

double L = double(12e-3); //Stage length of mine
/*
double DischargeElectrode2 = 0.0586;//from nozzle to skimmer
double Skimmer2FrontOfFirstCoil = 46.2e-3;
double CoilThickness = 9.6e-3;
double Skimmer2FirstStage = Skimmer2FrontOfFirstCoil + CoilThickness/2 - L;//from skimmer to the horizontal of the potential well
double Nozzle2FirstStage = DischargeElectrode2Skimmer + Skimmer2FirstStage;
*/

double CoilThickness = 9.6e-3;


//double Nozzle2Skimmer = 65.95e-3;
double Nozzle2Skimmer = 77.5e-3;
//65.95e-3 before installing adjustable nozzle on Nov 25, 2013
double SkimmerDiameter = 2e-3;
//double Nozzle2EntranceSlower = 157.23e-3;
//double Nozzle2EntranceSlower = 0.1508;
double Nozzle2EntranceSlower = 0.1508;//
//157.23e-3 before installing adjustable nozzle on Nov 25, 2013 
//not only Nozzle2Skimmer is changed, but whole decelerator moved downstream ~2mm

////////////parameter for matching simulation to experiment

double ShiftPosition = 0e-3;    //  Originally it was 0e-3  (24/July/2016)


double Nozzle2FirstStage = Nozzle2EntranceSlower + CoilThickness/2 - L/2 - ShiftPosition;


/*******Laser Detection Parameters*******/
/////////important for adjust the arrival time in simulation
//double DetectionPosition = 1.177;	//1.141, 1.138 used for trap timing simulation
double DetectionPosition = 1.177;		//Change dby Manish oin June 6th 2017
//after installing the adjustable nozzle, 1.141 for 80 stages


double IrisRadius = 0.0006; //0.0005
double ReighleyLength = 0.004; //0.002

double TransverseBound = 1.9e-3;	//slightly smaller than the inner radius of Zeeman Coil
 
//beam impossibly forms at the begining of nozzle opening time, i.e. time origin t=0, but is likely to form at the middle of nozzle opening 
double ValveOpeningTime = 60e-6;
double NO = ValveOpeningTime * pow(10,6);

double CoilOffestDelay = 0;//300e-6
double CD = CoilOffestDelay*pow(10,6);

double DischargeDelay = 60e-6;				//assume real radical is produced in 2nd pulse centered at 260us
double DD = DischargeDelay*pow(10,6);

double DischargeEnd = 160e-6;
double DE = DischargeEnd*pow(10,6);

double DischargePeak = 110e-6;		//110us initially
double DP = DischargePeak*pow(10,6);

//for new switching current profile with 80 coils on April 18,2013
//for SCR precharge time 20us, coil on time 80us from "coil global variables.vi", so IGBT on time 20+80=100us, SCR on time 80us
//delay from SCR on to current on is 13~16us
double OnDuration = 68e-6;

//double OnDuration = 80e-6;	//extend it to 80us to see how deceleration works
//several microsecond change will not affect the final velocity
double OffDelay = 7.5e-6;//3e-6, this will significantly affact the final velocity
double falling=OffDelay*pow(10,6);
//consider detailed timing taking account of IGBT and SCR
double DelayOfSwitchingOff = 3e-6;//time delay from falling edge of IGBT or SCR, i.e. the calculated switching-off time to real time of current switching-off
//this will not strongly affect the final velocity changing it from 2 to 6us
///////////////////////////////////////////

///////////////////////////////////////////

double dtSeqGen = 1e-6;
double DeltaT = 1e-6;
double DeltaTSeqGen = 1e-7;

double *MoleculeList=NULL;

/*
double t_start=2.6e-3;
double t_end=4.2e-3;

double t_start=4.18e-3;
double t_end=4.5e-3;

double t_start=4.45e-3;
double t_end=4.8e-3;

double t_start=4.94e-3;
double t_end=5.20e-3;

double t_start=5.50e-3;
double t_end=5.80e-3;

double t_start=5.84e-3;
double t_end=6.2e-3;

*/

double t_start=1.8e-3;
double t_end=8e-3;


int main(int argc, char* argv[])
{	

//////////////////////////////////////////////////////////////
	start_time = clock();
	
	srand( (unsigned)time( NULL ));


//write the magnitude of magnetic field into array

	int i,j,k;


	double TimeArray[100];
	double TimeOn[100];
	double VExit[100];

	for(i=0;i<100;i++)
	{
		TimeOn[i]=0;
		TimeArray[i]=0;
		VExit[i]=0;
	}

//	char* strTime= {"CH3_phi30.000000_v420.0mps_nozzle 60us_DP120us_DD 60us_DE170us_Falling  5us.txt"};

//	char* strTime= {"CH3_phi40.000000_v420.0mps_nozzle 60us_DP120us_DD 60us_DE170us_Falling  5us.txt"};
//	char* strTime= {"CH3_phi45.000000_v420.0mps_nozzle 60us_DP120us_DD 60us_DE170us_Falling  5us.txt"};
//	char* strTime= {"CH3_phi50.000000_v420.0mps_nozzle 60us_DP120us_DD 60us_DE170us_Falling  5us.txt"};
//	char* strTime= {"CH3_phi51.200000_v420.0mps_nozzle 60us_DP120us_DD 60us_DE170us_Falling  5us.txt"};

/*
	ReadTiming(strTime,TimeArray,TimeOn,VExit);

	for(i=0;i<NumberOfStages;i++)
	{
		TimeOn[i]+=CoilOffestDelay;	
		TimeArray[i]+=CoilOffestDelay;
	}

    cout<<"timing sequence has been loaded"<<endl;
*/

//initialization of the gradient of the field

	

	double *GradB = new double [ZNumber*RNumber*RNumber*3];

//loading gradient of the magnetic field along X direction
	char* DerivBx= {"dBbdX.txt"};

	ifstream fin3(DerivBx, ios::in);               
	if(fin3.fail()!=0)
	{
		cout<<"open data file failed"; 
		return 999 ;
	}

	for (i = 0; i < ZNumber; i++)
	{
		for (j = 0; j < RNumber; j++)
		{
			for (k = 0; k < RNumber; k++)
			{
				fin3 >> GradB[i*RNumber*RNumber*3+j*RNumber*3+k*3+0];
			}
		}
	}
	fin3.close();

//	cout<<"the gradient of the magnetic field along x direction has been successfully loaded into the array"<<endl;

//loading gradient of the magnetic field along Y direction
	char* DerivBy= {"dBbdY.txt"};

	ifstream fin4(DerivBy, ios::in);               
	if(fin4.fail()!=0)
	{
		cout<<"open data file failed"; 
		return 999 ;
	}

	for (i = 0; i < ZNumber; i++)
	{
		for (j = 0; j < RNumber; j++)
		{
			for (k = 0; k < RNumber; k++)
			{
				fin4 >> GradB[i*RNumber*RNumber*3+j*RNumber*3+k*3+1];
			}
		}
	}
	fin4.close();

//	cout<<"the gradient of the magnetic field along y direction has been successfully loaded into the array"<<endl;

//loading gradient of the magnetic field along Z direction

	char* DerivBz= {"dBbdZ.txt"};

	ifstream fin5(DerivBz, ios::in);               
	if(fin4.fail()!=0)
	{
		cout<<"open data file failed"; 
		return 999 ;
	}

	for (i = 0; i < ZNumber; i++)
	{
		for (j = 0; j < RNumber; j++)
		{
			for (k = 0; k < RNumber; k++)
			{
				fin5 >> GradB[i*RNumber*RNumber*3+j*RNumber*3+k*3+2];
			}
		}
	}
	fin5.close();

//	cout<<"the gradient of the magnetic field along z direction has been successfully loaded into the array"<<endl;

//loading fitting coefficient of the current profiles for 80 coils
//	char* fit= {"80CoilsCurrentProfileFit_June9thmod.txt"};
//	char* fit= {"CoilsCurrentProfileFit_Dec22.txt"};
//	char* fit= {"Try_CoilsCurrentProfileFit_Dec23.txt"};
//	char* fit= {"Dec24thFittingCurrentProfile.txt"};
//	char* fit= {"CurrentPolyFitting_Jan13.txt"};
	char* fit= {"Currentfitting18March2019.txt"};

	ifstream fin8(fit, ios::in);               
	if(fin8.fail()!=0)
	{
		cout<<"open data file failed"; 
		return 999 ;
	}

	double *fitting = new double [85*7];
	//80 coils totally, fitted with 4th order of polynomials

	for (i = 0; i < 85; i++)
	{
		for (j = 0; j < 7; j++)
		{
			fin8 >> fitting[i*7+j];
		}
	}
	fin8.close();



/////last for permanent magnets

	char c;
	int num=0;


	double *Bx_PM = new double [ZNumber0*YNumber0*XNumber0];
	num=0;
	ifstream is6;
	is6.open ("PMonlyBx.txt");
	while (is6.good())
	{
		c = is6.get();
		if (c=='\n')num++;
		is6 >> Bx_PM[num];	
	}
	is6.close();

	cout<<"Loading field of permanent magnets along x direction has been constructed!"<<endl;


	double *By_PM = new double [ZNumber0*YNumber0*XNumber0];
	num=0;
	ifstream is7;
	is7.open ("PMonlyBy.txt");
	while (is7.good())
	{
		c = is7.get();
		if (c=='\n')num++;
		is7 >> By_PM[num];	
	}
	is7.close();

	cout<<"Loading field of permanent magnets along y direction has been constructed!"<<endl;

	double *Bz_PM = new double [ZNumber0*YNumber0*XNumber0];
	num=0;
	ifstream is8;
	is8.open ("PMonlyBz.txt");
	while (is8.good())
	{
		c = is8.get();
		if (c=='\n')num++;
		is8 >> Bz_PM[num];	
	}
	is8.close();

	cout<<"Loading field of permanent magnets along z direction has been constructed!"<<endl;

//////////////////////





	double *Time;
	Time			=		new double[100000];
	int *ToFSignal;
	ToFSignal = new int[100000];

	for(i=0;i<100000;i++)
	{
		ToFSignal[i]=0;
	}

//	double TimeArray[800];
//	double VExit[800];
	char strTime[150];
//	sprintf(strTime,"TimeSequence_%2.1fT_phi%2.0f_v%3.1fmps_nozzle%3.0fus.txt",field,Phi,Vzcenter,NozzleOpenInUs);
	sprintf(strTime,"phi%2.2f_v%3.1fmps_nozzle%3.0fus_DD%3.0fus_DP%3.0fus.txt",Phi,Vsyn,NO,DD,DP);

//Generate the Timing Sequence
//	TimeArrayGenerator(Bn,DiffB,GradB,TimeArray,VExit,strTime);
	double t_arrive;
	NewFiniteSwitchingSequence(GradB,fitting,Bx_PM,By_PM,Bz_PM,TimeOn,TimeArray,VExit,&t_arrive,strTime);

//	traceing(Bn,DiffB,GradB);
 //   cout<<"time array has been generated"<<endl;
//	TraceWithSwitching(TimeArray,Bn,DiffB,GradB);	//trace 3 particles
//	trace(TimeArray,Bn,DiffB,GradB);	//trace single particle

//initialization of the array for the molecules when you run simulation for large number of molecules
/*
   if(MoleculeList!=NULL)delete MoleculeList;
	MoleculeList = new double [MoleculeNumber*9];

//assign position and velocity to molecules
	MoleculeEnsemble(MoleculeList);
   cout<<"moleculelist array has been initialized"<<endl;
*/

//	save_initial();

	char strTOF[50000];
	sprintf(strTOF,"Zeeman_TOF_In_B%2.1fT_v%dmps_Valve%3.0fus_Discharge%3.0fus.txt",field,Vsyn,NO,DP);

//	char strTOF2[50000];
//	sprintf(strTOF2,"Zeeman_TOF_Out_B%2.1fT_v%3.1fmps_Valve%3.0fus_ValveDela%3.0fus.txt",field,Vsyn,NO,ND);

	char strTOF3[50000];
	sprintf(strTOF3,"TOF_FreeFlight_v%3.0fmps_Valve%3.0fus_ValveDela%3.0fus.txt",vZcenter,NO,DP);
//	trajectory(Bn,DiffB,GradB);

/*
	double aa,bb;
	aa=DerivativeInterpolation(DiffB,55.5*1e-4,1,0,0);
	bb=DerivativeInterpolation(DiffB,1000*1e-4,1,0,0);
//	cout<<aa<<"		"<<bb<<endl;
	printf("%15.13e	%15.13e\n", aa, bb);
*/

//	FieldInterpolation(1,Bn,DiffB,GradB,-0.581e-3,0.239e-3,6.2e-3,1,0,0);

//Fly the Ensemble through the slower 
    New_SlowerFlight(TimeOn,TimeArray,Time,ToFSignal,GradB,fitting,strTOF);
	 cout<<"molecules has been slowed"<<endl;

//see the Time-of-flight spectrum for free molecules
//	FreeFly(Time,ToFSignal,strTOF3);
//	trace(TimeArray,Bn);


//delete allocated momory

	delete [] Bn;
	delete [] Shift;
	delete [] DiffB;
	delete [] GradB;
	delete[] Time;

	if(ToFSignal!=NULL)
	{
		delete ToFSignal;
		ToFSignal=NULL;
	}

	if(MoleculeList!=NULL)
	{
		delete MoleculeList;
		MoleculeList=NULL;
	}

	// end

	elapsed = clock()-start_time;

	elapsed_time = elapsed / ((double) CLOCKS_PER_SEC);
	printf("Time elapsed: %f\n",elapsed_time);

	return 0;
}

void MoleculeEnsemble(double *MoleculeList)
{
	//Declare variables to use later
	// Initializations of array
	long NumberAccepted;
	double Prob,PossibleDev;
	
	int i,j;

//	double DeltaZ = vZcenter*DischargeWidth/6;
/******Distribution Parameters *******/

  double Centers[6] = {Xcenter,Ycenter,Zcenter,vXcenter,vYcenter,vZcenter};
  double Deltas[6] = {DeltaX,DeltaY,DeltaZ,DeltavX,DeltavY,DeltavZ};

//initializing moleculeList

     for(i=0;i<6; i++)
	 {
	     NumberAccepted = 0;

		 while(NumberAccepted < MoleculeNumber)
		 {
//			 PossibleDev =  10*(double(rand())/RAND_MAX-0.5)*Deltas[i]; 
			 PossibleDev =  6*(double(rand())/RAND_MAX-0.5)*Deltas[i]; 
			 Prob = exp(-(PossibleDev/sqrt(2)/Deltas[i])*(PossibleDev/sqrt(2)/Deltas[i]));
			 if(double(rand())/RAND_MAX <= Prob)  
			 {
				 MoleculeList[NumberAccepted*9+i] = PossibleDev + Centers[i];
				 NumberAccepted += 1;
			 }

		 }
	 }

	double SumOfRatio[4]={0.33333,0.66666,0.83333,1.00000};//N(1)/N(0)=1/2
	NumberAccepted = 0;
	double my_tem;

	while(NumberAccepted < MoleculeNumber)
	{
		my_tem = double (rand())/RAND_MAX;
		MoleculeList[NumberAccepted*9+6] = -1;

		for(i=0;i<4;i++)
		{
			if(i==0)
			{
				if(my_tem >= 0 && my_tem <= SumOfRatio[i])
				{
					MoleculeList[NumberAccepted*9+6] = i;
				}
			}
			else
			{
				j=i-1;
				if(my_tem >= SumOfRatio[j] && my_tem <= SumOfRatio[i]) 
				{
					MoleculeList[NumberAccepted*9+6] = i;
				}
			}
			if(MoleculeList[NumberAccepted*9+6] != -1)
			{
				NumberAccepted += 1;
				break;
			}	
		}
	}

	NumberAccepted = 0;

	while(NumberAccepted < MoleculeNumber)
	{
		MoleculeList[NumberAccepted*9+7] = DischargePeak;	//time flag
		MoleculeList[NumberAccepted*9+8] = 0;	//position symbol
		NumberAccepted += 1;
	
	}

}


void save_initial()
{
	int i=0;
	char filename[30];
	sprintf(filename,"initialdata.txt");

	ofstream fout_initial(filename, ios::out);          //output initial molecule data
	if(fout_initial.fail()!=0)
	{
		cout<<"open failed"; 
		return;
	} 
	for(i=0;i<MoleculeNumber;i++) 
	{
		//counting numbers of molecules being loaded
	    fout_initial<<setw(14)<<MoleculeList[i*9+0]
			<<setw(14)<<MoleculeList[i*9+1]
			<<setw(14)<<MoleculeList[i*9+2]
			<<setw(14)<<MoleculeList[i*9+3]
			<<setw(14)<<MoleculeList[i*9+4]
			<<setw(14)<<MoleculeList[i*9+5]
			<<setw(14)<<MoleculeList[i*9+6]
			<<setw(14)<<MoleculeList[i*9+7]
			<<endl;  	
	}		
	fout_initial.close(); 
}


double PartialEnergyToB(double n)
{
	double uPrime;//for CH3
	//following parameters are fitted within 7 Tesla magnetic field
	if(n==0)//High-field-seeking state of |N=0,K=0,S=1/2,J=1/2>
	{uPrime=-0.46741*h*c*100;}
	if(n==1)//Low-field-seeking state of |N=0,K=0,S=1/2,J=1/2>
	{uPrime=0.46741*h*c*100;}
	if(n==2)//High-field-seeking state of |N=1,K=1,S=1/2,J=1/2>
	{uPrime=-0.4667*h*c*100;}
	if(n==3)//Low-field-seeking state of |N=1,K=1,S=1/2,J=3/2>
	{uPrime=0.46695*h*c*100;}
	return uPrime;
}


void MToAcceleration(double coeff,double *GradB,double x,double y,double z,double N,double *ax, double *ay, double *az)
{	
	int i,j,k,m;
	double dx,dy,dz,diff;
	double c[8],d[3];

	double xa,ya,xb,yb;

	i=(int)(floor)((x-RMin)/RGridSize);
	j=(int)(floor)((y-RMin)/RGridSize);
	k=(int)(floor)((z-ZMin)/ZGridSize);
	
	xa=RMin+i*RGridSize;
	ya=RMin+j*RGridSize;
	xb=RMin+(i+1)*RGridSize;
	yb=RMin+(j+1)*RGridSize;

//if (x,y) position of any one of 8 grids is outside the boundary	
	if(xa*xa+ya*ya>=TransverseBound*TransverseBound || xa*xa+yb*yb>=TransverseBound*TransverseBound || xb*xb+ya*ya>=TransverseBound*TransverseBound || xb*xb+yb*yb>=TransverseBound*TransverseBound)
	{
		*ax=0;
		*ay=0;
		*az=0;
	}

	else
	{
		if(i>=0 && i<RNumber-1 && j>=0 && j<RNumber-1 && k>=0 && k<ZNumber-1)
		{
			dx=x-(RMin+i*RGridSize);
			dy=y-(RMin+j*RGridSize);
			dz=z-(ZMin+k*ZGridSize);

			for(m=0;m<3;m++)
			{
				//linear interpolation of derivative of the magnetic field to get force
				c[0] = GradB[i*3+j*RNumber*3+k*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[1] = GradB[(i+1)*3+j*RNumber*3+k*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[2] = GradB[(i+1)*3+(j+1)*RNumber*3+k*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[3] = GradB[i*3+(j+1)*RNumber*3+k*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[4] = GradB[i*3+j*RNumber*3+(k+1)*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[5] = GradB[(i+1)*3+j*RNumber*3+(k+1)*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[6] = GradB[(i+1)*3+(j+1)*RNumber*3+(k+1)*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;
				c[7] = GradB[i*3+(j+1)*RNumber*3+(k+1)*RNumber*RNumber*3+m]*FieldScaleFactor*coeff;

				Interpolation(c,dx,dy,dz,&(d[m]));
				//d[] is interpolated magnetic field gradient of the specific point
			}

			diff=PartialEnergyToB(N);			

			*ax=(-1)*diff*d[0]/mass;
			*ay=(-1)*diff*d[1]/mass;
			*az=(-1)*diff*d[2]/mass;
		}
		else
		{
			*ax=0;
			*ay=0;
			*az=0;
		}
	}
}

int Interpolation(double *a,double dx,double dy,double dz, double *b)
{
	double ab,ef,gh,cd,aa,bb;

//	if(dx<0||dy<0||dz<0)return 0;
	if(dx<-pow(10,-10)||dy<-pow(10,-10)||dz<-pow(10,-10)) return 0;

	ab=a[0]+(a[1]-a[0])*dx/RGridSize;
	ef=a[4]+(a[5]-a[4])*dx/RGridSize;
	gh=a[7]+(a[6]-a[7])*dx/RGridSize;
	cd=a[3]+(a[2]-a[3])*dx/RGridSize;

	aa=ab+(cd-ab)*dy/RGridSize;
	bb=ef+(gh-ef)*dy/RGridSize;

	*b=aa+(bb-aa)*dz/ZGridSize;	//bn value

	return 1;
}

//coeff is a factor for scaling when current switching is on rising/falling edge
//n is a label for the state
void Runge_Kutta(double coeff,double *GradB,double x0,double y0,double z0,double vx0,double vy0,double vz0,double N,double *x,double *y,double *z,double *vx,double *vy,double *vz,double dt)
{
//define fi(t,vx,x,vy,y,vz,z)
//this method has higher precision than Euler(Newtonian)

	double k11,k12,k13,k14;
	double k21,k22,k23,k24;
	double k31,k32,k33,k34;
	double k41,k42,k43,k44;
	double k51,k52,k53,k54;
	double k61,k62,k63,k64;

	double ax,ay,az;
	double xx,yy,zz;
	
	if(coeff==0	|| x0*x0+y0*y0>=TransverseBound*TransverseBound)
	{
		*vx=vx0;
		*vy=vy0;
		*vz=vz0;
		*x=x0+vx0*dt;
		*y=y0+vy0*dt;
		*z=z0+vz0*dt;
	}
	else
	{
		MToAcceleration(coeff,GradB,x0,y0,z0,N,&ax,&ay,&az);
		k11=ax;//d(vx)/dt
		k21=vx0;//d(x)/dt
		k31=ay;//d(vy)/dt
		k41=vy0;//d(y)/dt
		k51=az;//d(vz)/dt
		k61=vz0;//d(z)/dt

		k22=vx0+dt/2*k11;	//medial vx
		xx=x0+dt/2*k21;		//medial x
		k42=vy0+dt/2*k31;	//medial vy
		yy=y0+dt/2*k41;		//medial y
		k62=vz0+dt/2*k51;	//medial vz
		zz=z0+dt/2*k61;		//medial z

		MToAcceleration(coeff,GradB,xx,yy,zz,N,&ax,&ay,&az);
		k12=ax;
		k32=ay;
		k52=az;

		k23=vx0+dt/2*k12;
		xx=x0+dt/2*k22;	
		k43=vy0+dt/2*k32;	
		yy=y0+dt/2*k42;	
		k63=vz0+dt/2*k52;	
		zz=z0+dt/2*k62;	

		MToAcceleration(coeff,GradB,xx,yy,zz,N,&ax,&ay,&az);
		k13=ax;
		k33=ay;
		k53=az;
			
		k24=vx0+dt*k13;	
		xx=x0+dt*k23;
		k44=vy0+dt*k33;	
		yy=y0+dt*k43;
		k64=vz0+dt*k53;	
		zz=z0+dt*k63;	

		MToAcceleration(coeff,GradB,xx,yy,zz,N,&ax,&ay,&az);
		k14=ax;
		k34=ay;
		k54=az;

		*vx=vx0+dt/6*(k11+2*k12+2*k13+k14);
		*x=x0+dt/6*(k21+2*k22+2*k23+k24);	
		*vy=vy0+dt/6*(k31+2*k32+2*k33+k34);
		*y=y0+dt/6*(k41+2*k42+2*k43+k44);
		*vz=vz0+dt/6*(k51+2*k52+2*k53+k54);
		*z=z0+dt/6*(k61+2*k62+2*k63+k64);

	}
}

//Defines Flight through the slower
void New_SlowerFlight(double TimeOn[],double TimeArray[],double Time[],int ToFSignal[],double *GradB,double *fitting,char *ptwriteZeemanToF)
{
	double x,y,z,vx,vy,vz,z0;
	double ttt=t_end + 0.5e-3;
	//	double SlowerEnd = Skimmer2FirstStage + L*double(NumberOfStages+1);
	int i,A;
	double coeff;
//	int n,nMin,nMax,temp_A,temp_B;
	int k=0;
	int kk=0;
	double step = 5e-6;
//	double random,random2,prob,Mb;
	double t5=TimeArray[5-1]+DelayOfSwitchingOff+OffDelay;
	double s5=Nozzle2EntranceSlower+(5-1)*L+CoilThickness/2;

	double t10=TimeArray[10-1]+DelayOfSwitchingOff+OffDelay;
	double s10=Nozzle2EntranceSlower+(10-1)*L+CoilThickness/2;
	
	double t20=TimeArray[20-1]+DelayOfSwitchingOff+OffDelay;
	double s20=Nozzle2EntranceSlower+(20-1)*L+CoilThickness/2;

		
	double t25=TimeArray[25-1]+DelayOfSwitchingOff+OffDelay;
	double s25=Nozzle2EntranceSlower+(25-1)*L+CoilThickness/2;
	
	double t30=TimeArray[30-1]+DelayOfSwitchingOff+OffDelay;
	double s30=Nozzle2EntranceSlower+(30-1)*L+CoilThickness/2;

	double t35=TimeArray[35-1]+DelayOfSwitchingOff+OffDelay;
	double s35=Nozzle2EntranceSlower+(35-1)*L+CoilThickness/2;

	double t40=TimeArray[40-1]+DelayOfSwitchingOff+OffDelay;
	double s40=Nozzle2EntranceSlower+(40-1)*L+CoilThickness/2;

	double t45=TimeArray[45-1]+DelayOfSwitchingOff+OffDelay;
	double s45=Nozzle2EntranceSlower+(45-1)*L+CoilThickness/2;
	
	double t50=TimeArray[50-1]+DelayOfSwitchingOff+OffDelay;
	double s50=Nozzle2EntranceSlower+(50-1)*L+CoilThickness/2;

	double t55=TimeArray[55-1]+DelayOfSwitchingOff+OffDelay;
	double s55=Nozzle2EntranceSlower+(55-1)*L+CoilThickness/2;

	double t60=TimeArray[60-1]+DelayOfSwitchingOff+OffDelay;
	double s60=Nozzle2EntranceSlower+(60-1)*L+CoilThickness/2;

	double t65=TimeArray[65-1]+DelayOfSwitchingOff+OffDelay;
	double s65=Nozzle2EntranceSlower+(65-1)*L+CoilThickness/2;
	
	double t70=TimeArray[70-1]+DelayOfSwitchingOff+OffDelay;
	double s70=Nozzle2EntranceSlower+(70-1)*L+CoilThickness/2;

	double t75=TimeArray[75-1]+DelayOfSwitchingOff+OffDelay;
	double s75=Nozzle2EntranceSlower+(75-1)*L+CoilThickness/2;

	double t80=TimeArray[80-1]+DelayOfSwitchingOff+OffDelay;
	double s80=Nozzle2EntranceSlower+(80-1)*L+CoilThickness/2;
	
	double t85=TimeArray[85-1]+DelayOfSwitchingOff+OffDelay;
	double s85=Nozzle2EntranceSlower+(85-1)*L+CoilThickness/2;

	ofstream fout_TOF(ptwriteZeemanToF, ios::out);          //open file for write
	if(fout_TOF.fail()!=0){cout<<"open starkTOF failed"; return;} 

	//Set the ToF time and initialize the signal array for this time

	//Increment for each slowing stage


	while(MoleculeList[1*9+7] < ttt)
	{
		Time[k] = MoleculeList[1*9+7];
		ToFSignal[k] = 0;  

		if(k%1000==0)	
		{
			cout<<"Now TOF is"<<"    "<<Time[k]<<endl;
		}		

//if molecules' time is before the switching sequence

		for(i = 0; i < MoleculeNumber; i++)
		{
			if(MoleculeList[i*9+8]>=0)
			{
				//searching for the molecule's position, which coil is this molecule inside 
				A = int(floor( (MoleculeList[i*9+2] - Nozzle2FirstStage)/L ));
				//B = int(ceil( (MoleculeList[i*9+2] - Nozzle2FirstStage)/L ));

				//first find out whether the molecule is inside the Zeeman decelerator or not

				if(A<0 || A>NumberOfStages-1 )
				{
					//molecule is outside the decelerator
					MoleculeList[i*9+0] += MoleculeList[i*9+3]*DeltaT;
					MoleculeList[i*9+1] += MoleculeList[i*9+4]*DeltaT;
					MoleculeList[i*9+2] += MoleculeList[i*9+5]*DeltaT;
				}
	//then finding out the current in which coil is turned on according to timing sequence
	//probably there is temporal overlap for current pulses in neighbourring coils
	//but due to the fact that only when the current on the coil A where the molecule lies is turned on,
	//will there be a force exerting on the molecule
				else
				{
					//molecule is inside the decelerator
					coeff=AppFittingSwitchingField(fitting,A,MoleculeList[i*9+7],TimeOn[A],TimeArray[A]);

					//the following program will be the same whether coeff==0 or not
					//since Runge_Kutta function already consider both cases
					z0 = MoleculeList[i*9+2] - A*L - Nozzle2FirstStage;
					Runge_Kutta(coeff,GradB,MoleculeList[i*9+0],MoleculeList[i*9+1],z0,MoleculeList[i*9+3],MoleculeList[i*9+4],MoleculeList[i*9+5],MoleculeList[i*9+6],&x,&y,&z,&vx,&vy,&vz,DeltaT);
					MoleculeList[i*9+0] = x;
					MoleculeList[i*9+1] = y;
					MoleculeList[i*9+2] = z + A*L + Nozzle2FirstStage;
					MoleculeList[i*9+3] = vx;
					MoleculeList[i*9+4] = vy;
					MoleculeList[i*9+5] = vz;


	//and if coil B or coil C is on at this time, they may have small interaction with the molecule 
	//because of small penetrating field produced by the coil at the molecule's position
	//but here we just ignore it since the calculated magnetic field is only for 1*L instead of 2*L	
				}
			}
				
			MoleculeList[i*9+7] += DeltaT;
			
			if(sqrt(MoleculeList[i*9]*MoleculeList[i*9]+MoleculeList[i*9+1]*MoleculeList[i*9+1])>=TransverseBound)
			{
				MoleculeList[i*9+8]=-1;
			}
			if(MoleculeList[i*9+8]>=0 &&  (MoleculeList[i*9+2] - DetectionPosition)*(MoleculeList[i*9+2] - DetectionPosition)+MoleculeList[i*9]*MoleculeList[i*9] <= IrisRadius*IrisRadius && fabs(MoleculeList[i*9+1]) < ReighleyLength)   
			{ 
				ToFSignal[k]+=1;

			}
		} //end of for loop of molecules

//for obtaining phase-space plot
	/*	
		if(fabs(Time[k]-t5)<DeltaT)
		{
			end_on(t5,s5);
		}

		if(fabs(Time[k]-t10)<DeltaT)
		{
			end_on(t10,s10);
		}
		if(fabs(Time[k]-t20)<DeltaT)
		{
			end_on(t20,s20);
		}
		if(fabs(Time[k]-t25)<DeltaT)
		{
			end_on(t25,s25);
		}
		if(fabs(Time[k]-t30)<DeltaT)
		{
			end_on(t30,s30);
		}
		if(fabs(Time[k]-t35)<DeltaT)
		{
			end_on(t35,s35);
		}

		if(fabs(Time[k]-t40)<DeltaT)
		{
			end_on(t40,s40);
		}
		if(fabs(Time[k]-t45)<DeltaT)
		{
			end_on(t45,s45);
		}
		if(fabs(Time[k]-t50)<DeltaT)
		{
			end_on(t50,s50);
		}
		if(fabs(Time[k]-t55)<DeltaT)
		{
			end_on(t55,s55);
		}
		if(fabs(Time[k]-t60)<DeltaT)
		{
			end_on(t60,s60);
		}
		if(fabs(Time[k]-t65)<DeltaT)
		{
			end_on(t65,s65);
		}
		if(fabs(Time[k]-t70)<DeltaT)
		{
			end_on(t70,s70);
		}
		if(fabs(Time[k]-t75)<DeltaT)
		{
			end_on(t75,s75);
		}
		if(fabs(Time[k]-t80)<DeltaT)
		{
			end_on(t80,s80);
		}

		if(fabs(Time[k]-t85)<DeltaT)
		{
			end_on(t85,s85);
		}
*/
//		fout_TOF<<1000*Time[k]<<"     "<<ToFSignal[k]<<endl; 

		//Here 5 is the scan step for free flight
		if( MoleculeList[1*9+7]>=t_start && MoleculeList[1*9+7]<=ttt && int(MoleculeList[1*9+7]*1e6) % 5 == 0 )
		{
			fout_TOF<<setw(14)<<1000*Time[k]
					<<setw(14)<<ToFSignal[k]
					<<endl;      //output TOF signal to file 

		}

		k += 1;
	}

	fout_TOF.close();
	cout<<"Output to file stark_TOF successfully!"<<endl;

}

void FreeFly(double Time[],int ToFSignal[],char *ptwriteZeemanToF)
{
	int i=0,k=0;
	int KMax = 1000;
	double step = 20e-6;

	//open file for TOF signal output
	ofstream fout_TOF(ptwriteZeemanToF, ios::out);          //open file for write
	if(fout_TOF.fail()!=0){cout<<"open starkTOF failed"; return;} 

	while(MoleculeList[1*9+7]<= t_end)
	{			
		if(k%500==0){cout<<"Now k is"<<"    "<<k<<endl;}

		Time[k] = MoleculeList[1*9+7];
		ToFSignal[k] = 0;   //counter the number of molecules reaching the Detection position after each timestep, or say each DeltaT

		for(int i = 0; i < MoleculeNumber; i++)
		{
			MoleculeList[i*9+0] += MoleculeList[i*9+3]*DeltaT;
			MoleculeList[i*9+1] += MoleculeList[i*9+4]*DeltaT;
			MoleculeList[i*9+2] += MoleculeList[i*9+5]*DeltaT;
			MoleculeList[i*9+7] += DeltaT;
				
			if(MoleculeList[i*9+2] < Nozzle2Skimmer && sqrt(MoleculeList[i*9]*MoleculeList[i*9]+MoleculeList[i*9+1]*MoleculeList[i*9+1])>=SkimmerDiameter/2)
			{
				MoleculeList[i*9+8]=-1;
			}
								
			if(MoleculeList[i*9+2] < Nozzle2EntranceSlower && sqrt(MoleculeList[i*9]*MoleculeList[i*9]+MoleculeList[i*9+1]*MoleculeList[i*9+1])>=TransverseBound)
			{
				MoleculeList[i*9+8]=-1;
			}

			if(MoleculeList[i*9+8]>=0 &&  (MoleculeList[i*9+2] - DetectionPosition)*(MoleculeList[i*9+2] - DetectionPosition)+MoleculeList[i*9]*MoleculeList[i*9] <= IrisRadius*IrisRadius && fabs(MoleculeList[i*9+1]) < ReighleyLength)   
			{ 
				ToFSignal[k]+=1;

			}

		}
//		fout_TOF<<1000*Time[k]<<"     "<<ToFSignal[k]<<endl;  

		if( MoleculeList[1*9+7]>=t_start && MoleculeList[1*9+7]<=t_end && int(MoleculeList[1*9+7]*1e6) % 20 == 0 )
		{
			fout_TOF<<setw(14)<<1000*Time[k]
					<<setw(14)<<ToFSignal[k]
					<<endl;      //output TOF signal to file 


		}

//		if(fabs(Time[k]-t1*1e-3)<0.5*DeltaT){save_TOF(t1);}
//		if(fabs(Time[k]-t2*1e-3)<0.5*DeltaT){save_TOF(t2);}

/*
		if(fabs(Time[k]-ts)<=DeltaT/2)
		{
			end_off(ts);      //output TOF signal to file 
		}
*/
		k += 1;
	}

	fout_TOF.close();

	cout<<"free flight end!"<<endl;
}


void end_off(double t)
{
	int i;
	double tt=t*1e3;

	char filename[30];
	sprintf(filename,"t=%1.3fOFF.txt",tt);

	ofstream fout_end(filename, ios::out);          //output initial molecule data
	if(fout_end.fail()!=0)
	{
		cout<<"open failed"; 
		return;
	} 

	for(i=0;i<MoleculeNumber;i++) 
	{
			 fout_end<<setw(14)<<MoleculeList[i*9+0]
						<<setw(14)<<MoleculeList[i*9+1]
						<<setw(14)<<MoleculeList[i*9+2]
						<<setw(14)<<MoleculeList[i*9+3]
						<<setw(14)<<MoleculeList[i*9+4]
						<<setw(14)<<MoleculeList[i*9+5]
						<<endl; 
	}
	
	fout_end.close();

}

//for recording phase space of molecules
void end_on(double t,double z)
{
	int i;
	double tt=t*1e3;

	char filename[30];
	sprintf(filename,"t=%1.3fON.txt",tt);

	ofstream fout_end(filename, ios::out);          //output initial molecule data
	if(fout_end.fail()!=0)
	{
		cout<<"open failed"; 
		return;
	} 

	for(i=0;i<MoleculeNumber;i++) 
	{
		if(fabs(MoleculeList[i*9+2]-z)<=L && fabs(MoleculeList[i*9+0])<2*TransverseBound && fabs(MoleculeList[i*9+1])<2*TransverseBound)
		{
			 fout_end<<setw(14)<<MoleculeList[i*9+0]
						<<setw(14)<<MoleculeList[i*9+1]
						<<setw(14)<<MoleculeList[i*9+2]
						<<setw(14)<<MoleculeList[i*9+3]
						<<setw(14)<<MoleculeList[i*9+4]
						<<setw(14)<<MoleculeList[i*9+5]
						<<setw(14)<<MoleculeList[i*9+6]
						<<endl; 
		}
	}
	
	fout_end.close();

}

double SwitchingField(double t,double t0)
{
	double ff;

//consider both rising time and falling timing

//new fitted function for 80 coils,300A on April 19,2013
/*
	double a=1.9585;
	double b=0.40348;
	double c=1.6101e-5;
*/

//new fitted function for 80 coils,600A on April 23,2013

	double a=2.097;
	double b=0.43794;
	double c=1.8681e-5;

	double t1=t0+DelayOfSwitchingOff;	//real switching-off time
	double t2=t1-OnDuration;	//real switching-on time
	double dt=t-t2;	//the time from current switching on

	double gg=a+b*log10(OnDuration+c);

	if(dt>0 && dt<OnDuration)
	{
		ff=a+b*log10(dt+c);
		ff/=gg;
	}
	if(dt>OnDuration && dt<OnDuration+OffDelay)
	{
		ff=1-(dt-OnDuration)/OffDelay;
	}
	if(dt<=0 || dt>=OnDuration+OffDelay)
	{
		ff=0;
	}

	return ff;

}

void FieldInterpolation(double coeff,double *GradB,double x,double y,double z)
{	
	double d[3];
	FILE * pFile;
	pFile = fopen ("interpolation.txt","w");

	int i,j,k,m;
	double dx,dy,dz;
	double c[8];

	double xa,ya,xb,yb;

	i=(int)(floor)((x-RMin)/RGridSize);
	j=(int)(floor)((y-RMin)/RGridSize);
	k=(int)(floor)((z-ZMin)/ZGridSize);
	
	xa=RMin+i*RGridSize;
	ya=RMin+j*RGridSize;
	xb=RMin+(i+1)*RGridSize;
	yb=RMin+(j+1)*RGridSize;

	dx=x-(RMin+i*RGridSize);
	dy=y-(RMin+j*RGridSize);
	dz=z-(ZMin+k*ZGridSize);

	for(m=0;m<3;m++)
	{
		//linear interpolation of derivative of the magnetic field to get force
		c[0] = GradB[i*3+j*RNumber*3+k*RNumber*RNumber*3+m];
		c[1] = GradB[(i+1)*3+j*RNumber*3+k*RNumber*RNumber*3+m];
		c[2] = GradB[(i+1)*3+(j+1)*RNumber*3+k*RNumber*RNumber*3+m];
		c[3] = GradB[i*3+(j+1)*RNumber*3+k*RNumber*RNumber*3+m];
		c[4] = GradB[i*3+j*RNumber*3+(k+1)*RNumber*RNumber*3+m];
		c[5] = GradB[(i+1)*3+j*RNumber*3+(k+1)*RNumber*RNumber*3+m];
		c[6] = GradB[(i+1)*3+(j+1)*RNumber*3+(k+1)*RNumber*RNumber*3+m];
		c[7] = GradB[i*3+(j+1)*RNumber*3+(k+1)*RNumber*RNumber*3+m];

		Interpolation(c,dx,dy,dz,&(d[m]));
		//d[] is interpolated magnetic field gradient of the specific point
	}

	fprintf (pFile, "%15.13e\t%15.13e\t%15.13e\n",d[0],d[1],d[2]);
	fclose (pFile);

}



double NewFittingSwitchingField(double *fitting,int i,double t,double ton)
{
	double ff;

	double k0=fitting[i*7+0];	//maximal pulse width ever fitted for this coil current profile, in units of us
	double k1=fitting[i*7+1];
	double k2=fitting[i*7+2];
	double k3=fitting[i*7+3];
	double k4=fitting[i*7+4];
	double k5=fitting[i*7+5];
	double k6=fitting[i*7+6];

	double dt=(t-ton)*1e6;	//real on time
	double t_error=5;

	double gg=k1+k2*dt+k3*dt*dt+k4*dt*dt*dt+k5*dt*dt*dt*dt+k6*dt*dt*dt*dt*dt;

	if(dt>0 && dt<k0+t_error)
	{
		ff=gg/600;
	}
	if(dt>=k0+t_error)
	{
		ff=0;
		cout<<"current pulse width exceeds the limitation because of fitting"<<endl;
	}

	return ff;

}


double AppFittingSwitchingField(double *fitting,int i,double t,double ton,double toff)
{
	double ff;

	double k0=fitting[i*7+0];
	double k1=fitting[i*7+1];
	double k2=fitting[i*7+2];
	double k3=fitting[i*7+3];
	double k4=fitting[i*7+4];
	double k5=fitting[i*7+5];
	double k6=fitting[i*7+6];

	double dt=(t-ton)*1e6;	//real on time
	double coil_OnTime = (toff-ton)*1e6;
//	double max_t=110e-6;	//on-time when current reaches the maximum ~630A

	double gg;
	double ga;

	if(dt<=0 || dt>=coil_OnTime+OffDelay)	//if(t<=ton+DelayOfSwitchingOff && t>=toff+DelayOfSwitchingOff+OffDelay)
	{
		ff=0;
	}

	else
	{
		if(dt>0 && dt<coil_OnTime)		//if(t>ton+DelayOfSwitchingOff && t<toff+DelayOfSwitchingOff)
		{
			gg=k1+k2*dt+k3*dt*dt+k4*dt*dt*dt+k5*dt*dt*dt*dt+k6*dt*dt*dt*dt*dt;
			ff=gg/600;
		}

		if(dt>=coil_OnTime && dt<coil_OnTime+OffDelay*1e6)		//if(t>=toff+DelayOfSwitchingOff && t<toff+DelayOfSwitchingOff+OffDelay)
		{
			ga=k1+k2*coil_OnTime+k3*coil_OnTime*coil_OnTime+k4*coil_OnTime*coil_OnTime*coil_OnTime+k5*coil_OnTime*coil_OnTime*coil_OnTime*coil_OnTime+k6*coil_OnTime*coil_OnTime*coil_OnTime*coil_OnTime*coil_OnTime;
			ff=(ga/600)*(toff+OffDelay-t)/OffDelay;
		}
	}

	return ff;
}


/*
//if the On_time is too long such that the current starts to decrease, then exponential fitting is not good
double NewFittingSwitchingField(double *fitting,int i,double t,double ton)
{
	double ff;

	double k0=fitting[i*3+0];
	double k1=fitting[i*3+1];
	double k2=fitting[i*3+2];

	double dt=(t-ton)*1e6;	//real on time

	double gg=k0+k1*exp(-k2*dt);	

	if(dt>0 && dt<200)
	{
		ff=gg/600;
	}

	return ff;

}


double AppFittingSwitchingField(double *fitting,int i,double t,double ton,double toff)
{
	double ff;

	double k0=fitting[i*3+0];
	double k1=fitting[i*3+1];
	double k2=fitting[i*3+2];

	double dt=(t-ton)*1e6;	//real on time
	double coil_OnTime = (toff-ton)*1e6;
//	double max_t=110e-6;	//on-time when current reaches the maximum ~630A

	double gg;
	double ga;

	if(dt<=0 || dt>=coil_OnTime+OffDelay)	//if(t<=ton+DelayOfSwitchingOff && t>=toff+DelayOfSwitchingOff+OffDelay)
	{
		ff=0;
	}

	else
	{
		if(dt>0 && dt<coil_OnTime)		//if(t>ton+DelayOfSwitchingOff && t<toff+DelayOfSwitchingOff)
		{
			gg=k0+k1*exp(-k2*dt);
			ff=gg/600;
		}

		if(dt>=coil_OnTime && dt<coil_OnTime+OffDelay*1e6)		//if(t>=toff+DelayOfSwitchingOff && t<toff+DelayOfSwitchingOff+OffDelay)
		{
			ga=k0+k1*exp(-k2*coil_OnTime);
			ff=(ga/600)*(toff+OffDelay-t)/OffDelay;
		}
	}

	return ff;
}
*/


double ReadTiming(char* PRead,double a[],double b[],double c[])
{
	int counter;
	ifstream fin(PRead, ios::in);               
	if(fin.fail()!=0)
	{
		cout<<"open data file failed"; 
		return 999 ;
	}
	counter=0;
	while(counter<NumberOfStages)
	{
		fin>>a[counter]>>b[counter]>>c[counter];
		counter++;
	}
	 fin.close();
	 return 0 ;

}


void NewFiniteSwitchingSequence(double *GradB,double *fitting,double *Bx_PM,double *By_PM,double *Bz_PM,double TimeOn[],double TimeArray[],double VExit[],double *ta,char* ptwriteTimeSequence)
{
	//Phi is in degrees.

	double x0 = 0;
	double y0 = 0;
	double z0 = Zcenter; //Initialize the molecule at the end of the hexapole
	double vx0 = 0;
	double vy0 = 0;
	double vz0 = Vsyn; //Initialize the molecule with vZcenter as the velocity
	double x,y,z,vx,vy,vz;
	double t = DischargePeak;
	double ZStop,ZLocal;
	int k=0;
	int k_temp=0;

	double t_ahead = 85e-6;

	double t_arrival;

//for the selected state
	int m=1;				//1 for Low-field-seeking state, 0 for high-field-seeking state
///////////////////

	double scaling=0;
/*
	char filename[30];
	sprintf(filename,"coeff.txt");

	ofstream fout_coeff(filename, ios::out);          //output initial molecule data
	if(fout_coeff.fail()!=0)
	{
		cout<<"open failed"; 
		return;
	} 
*/


//   FILE * pFile;
//   pFile = fopen ("coeff.txt","w");

	int A,B;
//	int ClosestStage,A,B,StageUse;

	ofstream fout(ptwriteTimeSequence, ios::out);          //open file for write
    if(fout.fail()!=0){cout<<"open failed"; return;} 

	for(int j = 0; j < NumberOfStages; j++)                   
	{
		ZStop = Nozzle2FirstStage + (Phi*L/180) + j*L; 

		//Fly molecules until fields are switched
		while( (z0 < ZStop) && (vz0 > 0) )
		{

		//This section finds the closest charged coil to the molecules, which depends on how coils are connected and switched on/off
		//if coils are connected to each other every N coils, this section is meaningful
		//if coils are independent, not connected to other coils every N coils, they are switched independently, then this section has no meaning
			A = int(floor( (z0-Nozzle2FirstStage)/L ));  
		//	B = int(ceil( z0/L ));   
 			B = A+1; 
					
			//If Molecules are not yet in the slower field or past the slower let them free fly
			//for the field array loaded, it ranges from 0 to 2L (-Pi/2 ~ Pi/2)
			if(B <= j)		//same as if(A<j)
			{
				z0 += vz0*DeltaTSeqGen;
				t += DeltaTSeqGen;
				k+=1;
				scaling=0;
				k_temp = k;
		//		fout_coeff <<t<<"   "<<scaling<<"	"<<z0<<"	"<<vz0<< endl;
		//		fprintf (pFile,"%10.6e	%10.6f\n", t, vz0);

			}
			else
			{
              //Now determine the local (in terms of the force array) coordinates
				//Find the DeltaV kick and propagate the molecule
				ZLocal = z0 - Nozzle2FirstStage - j*L;
				//if(fabs(ZLocal)<=L*1e-1)	TimeOn[j]=t-t_ahead;

				///////////////////////preset of Switching_On Time for each coil
				if(k==k_temp)	
				{
					t_arrival=t;
					TimeOn[j]=t_arrival-t_ahead;

					if(j>0 && TimeOn[j]>TimeArray[j-1])					//In order to temporally overlap neighbouring coils
					{
						TimeOn[j]=TimeArray[j-1];
						cout<<"For Stage	"<<j+1<<"  "<<"Current is Switched on In Advance to make Neighbouring Coils Temporally Overlap"<<endl;
					}

					///////// Initial and Rough Judging whether Preset Value of TimeOn[j] is good or bad, and reset it if it is bad

					if(Phi*L/180/vz0+(t_arrival-TimeOn[j])>=fitting[j*7]*1e-6)
					{
						cout<<"Switching_On Current early to Temporally Overlap Neighbouring Coils is NOT Suitable and Practical"<<endl;
						if(Phi*L/180/vz0<=fitting[j*7]*1e-6)
						{
							TimeOn[j]=t_arrival+Phi*L/180/vz0-fitting[j*7]*1e-6;
						}
						else
						{
							TimeOn[j]=t_arrival;
						}
					}
				}

				scaling=NewFittingSwitchingField(fitting,j,t,TimeOn[j]);
				Runge_Kutta(scaling,GradB,x0,y0,ZLocal,vx0,vy0,vz0,m,&x,&y,&z,&vx,&vy,&vz,dtSeqGen);

				x0=x;
				y0=y;
				z0=z+Nozzle2FirstStage+j*L;
				vx0=vx;
				vy0=vy;
				vz0=vz;
				t += dtSeqGen;
				k+=1;


			}
		}
		
		TimeArray[j] = t; //Get Time for field switching

		cout<<"MaxTime for Coil "<<"	"<<j+1<<"  "<<fitting[j*7]<<endl;
		cout<<"OnTime for Coil "<<"	"<<j+1<<"  "<<(TimeArray[j]-TimeOn[j])*1e6<<endl;
		cout<<"Final velocity "<<"	"<<vz0<<endl;


		//coil OnTime should not exceed the MaxOnTime used for current profile fitting,otherwise go back to the past when the target molecule enters this coil
		while((TimeArray[j]-TimeOn[j])*1e6>fitting[j*7])
		{
			t=t_arrival;					//the time when synchronous molecule reach the coil
			TimeOn[j]+=5e-6;				//increment the TimeOn by 1us

			///////////////////Switch on current before the synchronous molecule reach the coil
	
			x0=0;
			y0=0;
			z0=Nozzle2FirstStage+j*L;
			vx0=0;
			vy0=0;
			vz0=VExit[j-1];


			if(TimeOn[j]<=t_arrival)	//Switching-On Time is still ahead of Arrival Time after increment, then keep adjusted TimeOn[j]
			{
				while(z0-Nozzle2FirstStage-j*L < Phi*L/180)				//Switch off when syn molecule reach the position corresponding to nominated phase angle
				{
					ZLocal=z0-Nozzle2FirstStage-j*L;
					scaling=NewFittingSwitchingField(fitting,j,t,TimeOn[j]);
					Runge_Kutta(scaling,GradB,x0,y0,ZLocal,vx0,vy0,vz0,m,&x,&y,&z,&vx,&vy,&vz,dtSeqGen);

					x0=x;
					y0=y;
					z0=z+Nozzle2FirstStage+j*L;
					vx0=vx;
					vy0=vy;
					vz0=vz;
					t += dtSeqGen;
					k+=1;
				}
				TimeArray[j] = t;
			}
			else						//Switching-On Time now fall behind Arrival Time after increment, then adjust TimeOn[j] to be the Arrival Time
			{
				TimeOn[j]=t_arrival;				//increment the TimeOn by 1us
				cout<<"Switch on current when the synchronous molecule reach the coil"<<endl;


				//////////////////////Switch off either when syn molecule reach the position corresponding to nominated phase angle, or when OnTime reaches Max Duration of current pulse
				while(z0-Nozzle2FirstStage-j*L < Phi*L/180	|| t-TimeOn[j]<=fitting[j*7]*1e-6)
				{
					ZLocal=z0-Nozzle2FirstStage-j*L;
					scaling=NewFittingSwitchingField(fitting,j,t,TimeOn[j]);
					Runge_Kutta(scaling,GradB,x0,y0,ZLocal,vx0,vy0,vz0,m,&x,&y,&z,&vx,&vy,&vz,dtSeqGen);

					x0=x;
					y0=y;
					z0=z+Nozzle2FirstStage+j*L;
					vx0=vx;
					vy0=vy;
					vz0=vz;
					t += dtSeqGen;
					k+=1;
				}
				TimeArray[j] = t;
			}
		}

		cout<<"OnTime for Coil "<<"	"<<j+1<<"  "<<(TimeArray[j]-TimeOn[j])*1e6<<endl;
		cout<<"Final velocity "<<"	"<<vz0<<endl;
			cout<<"Final position "<<"	"<<z0<<endl;

		while(t>=TimeArray[j] && t<=TimeArray[j]+OffDelay && z0<Nozzle2FirstStage+(j+1)*L)
		{
				ZLocal = z0 - Nozzle2FirstStage - j*L;
				scaling=NewFittingSwitchingField(fitting,j,TimeArray[j],TimeOn[j])*(OffDelay-t+TimeArray[j])/OffDelay;
				Runge_Kutta(scaling,GradB,x0,y0,ZLocal,vx0,vy0,vz0,m,&x,&y,&z,&vx,&vy,&vz,dtSeqGen);

				x0=x;
				y0=y;
				z0=z+Nozzle2FirstStage+j*L;
				vx0=vx;
				vy0=vy;
				vz0=vz;
				t += dtSeqGen;
				k+=1;
			//	fout_coeff <<t<<"   "<<scaling<<"	"<<z0<<"	"<<vz0<< endl;
			//	fprintf (pFile,"%10.6e	%10.6f\n", t, vz0);

		}

		VExit[j] = vz0;

	}

	for(j=0;j<NumberOfStages;j++)
	{
		fout<<TimeArray[j]-DelayOfSwitchingOff<<"   "<<TimeOn[j]-DelayOfSwitchingOff<<"   "<<VExit[j]<<endl;	
	}

//	fout_coeff.close();
	fout.close(); 
	cout<<"Output to file TimeSequence successfully!"<<endl;
	cout<<"After decelerator the velocity is"<<"		"<<vz0<<endl;
////////////////////////

//////////////////////////After switching off last coil of decelerator, trace molecule trapped by permanent magnets only

	while(t>TimeArray[84]+OffDelay && z0-DetectionPosition<=ZMin0)
	{
		x0 += vx0*DeltaTSeqGen;
		y0 += vy0*DeltaTSeqGen;
		z0 += vz0*DeltaTSeqGen;
		t += DeltaTSeqGen;
	}

	while(t>TimeArray[84]+OffDelay && z0-DetectionPosition>ZMin0 && z0-DetectionPosition<=0)
	{
		ZLocal = z0-DetectionPosition;

		Trap_Runge_Kutta(Bx_PM,By_PM,Bz_PM,x0,y0,ZLocal,vx0,vy0,vz0,m,&x,&y,&z,&vx,&vy,&vz,dtSeqGen);

		x0 = x;
		y0 = y;
		z0 = z + DetectionPosition;
		vx0 = vx;
		vy0 = vy;
		vz0 = vz;
		t += dtSeqGen;
	}

	*ta=t;
	cout<<"the final velocity is"<<"		"<<vz0<<endl;
	cout<<"when does the synchronous molecule arrive the detector"<<"		"<<*ta<<endl;
}



void TrapAcceleration(double *Btx,double *Bty,double *Btz,double x,double y,double z,double m,double *ax, double *ay, double *az)
{	
	int i,j,k;
	double bx,by,bz,dx,dy,dz,db2dx,db2dy,db2dz;
	double diff;
	double a[8],b[8],c[8];

	i=(int)(floor)((x-XMin0)/GridSize0);
	j=(int)(floor)((y-YMin0)/GridSize0);
	k=(int)(floor)((z-ZMin0)/GridSize0);

	if(i>=0 && i<XNumber0-1 && j>=0 && j<YNumber0-1 && k>=0 && k<ZNumber0-1)
	{
		dx=x-(XMin0+i*GridSize0);
		dy=y-(YMin0+j*GridSize0);
		dz=z-(ZMin0+k*GridSize0);

		a[0]=fabs(Btx[k*YNumber0*XNumber0+j*XNumber0+i]);
		a[1]=fabs(Btx[k*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		a[2]=fabs(Btx[k*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		a[3]=fabs(Btx[k*YNumber0*XNumber0+(j+1)*XNumber0+i]);
		a[4]=fabs(Btx[(k+1)*YNumber0*XNumber0+j*XNumber0+i]);
		a[5]=fabs(Btx[(k+1)*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		a[6]=fabs(Btx[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		a[7]=fabs(Btx[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+i]);

		InterpolationX_trap(a,dx,dy,dz,&bx,&db2dx);

		b[0]=fabs(Bty[k*YNumber0*XNumber0+j*XNumber0+i]);
		b[1]=fabs(Bty[k*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		b[2]=fabs(Bty[k*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		b[3]=fabs(Bty[k*YNumber0*XNumber0+(j+1)*XNumber0+i]);
		b[4]=fabs(Bty[(k+1)*YNumber0*XNumber0+j*XNumber0+i]);
		b[5]=fabs(Bty[(k+1)*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		b[6]=fabs(Bty[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		b[7]=fabs(Bty[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+i]);

		InterpolationY_trap(b,dx,dy,dz,&by,&db2dy);

		c[0]=fabs(Btz[k*YNumber0*XNumber0+j*XNumber0+i]);
		c[1]=fabs(Btz[k*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		c[2]=fabs(Btz[k*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		c[3]=fabs(Btz[k*YNumber0*XNumber0+(j+1)*XNumber0+i]);
		c[4]=fabs(Btz[(k+1)*YNumber0*XNumber0+j*XNumber0+i]);
		c[5]=fabs(Btz[(k+1)*YNumber0*XNumber0+j*XNumber0+(i+1)]);
		c[6]=fabs(Btz[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+(i+1)]);
		c[7]=fabs(Btz[(k+1)*YNumber0*XNumber0+(j+1)*XNumber0+i]);

		InterpolationZ_trap(c,dx,dy,dz,&bz,&db2dz);
		//d is interpolated magnetic field strength of the specific point

		diff=PartialEnergyToB(m);			

		*ax=(-1)*diff*db2dx/mass;
		*ay=(-1)*diff*db2dy/mass;
		*az=(-1)*diff*db2dz/mass;

	}

	else
	{
		*ax=0;
		*ay=0;
		*az=0;
	}

}



int InterpolationZ_trap(double *a,double dx,double dy,double dz, double *b,double *c)
{
	double ab,ef,gh,cd,aa,bb;

//	if(dx<0||dy<0||dz<0)return 0;
	if(dx<-pow(10,-10)||dy<-pow(10,-10)||dz<-pow(10,-10)) return 0;

	ab=a[0]+(a[1]-a[0])*dx/GridSize0;
	ef=a[4]+(a[5]-a[4])*dx/GridSize0;
	gh=a[7]+(a[6]-a[7])*dx/GridSize0;
	cd=a[3]+(a[2]-a[3])*dx/GridSize0;

	aa=ab+(cd-ab)*dy/GridSize0;
	bb=ef+(gh-ef)*dy/GridSize0;

	*b=aa+(bb-aa)*dz/GridSize0;	//bn value
	*c=(bb-aa)/GridSize0;		//dbn2dz

	return 1;
}

int InterpolationY_trap(double *a,double dx,double dy,double dz, double *b,double *c)
{
	double ab,ef,gh,cd,cc,dd;

//	if(dx<0||dy<0||dz<0)return 0;
	if(dx<-pow(10,-10)||dy<-pow(10,-10)||dz<-pow(10,-10)) return 0;

	ab=a[0]+(a[1]-a[0])*dx/GridSize0;
	ef=a[4]+(a[5]-a[4])*dx/GridSize0;
	gh=a[7]+(a[6]-a[7])*dx/GridSize0;
	cd=a[3]+(a[2]-a[3])*dx/GridSize0;

	cc=ab+(ef-ab)*dz/GridSize0;
	dd=cd+(gh-cd)*dz/GridSize0;

	*b=cc+(dd-cc)*dy/GridSize0;	//bn value
	*c=(dd-cc)/GridSize0;		//dbn2dy

	return 1;
}

int InterpolationX_trap(double *a,double dx,double dy,double dz, double *b,double *c)
{
	double ae,dh,bf,cg,ee,ff;

//	if(dx<0||dy<0||dz<0)return 0;
	if(dx<-pow(10,-10)||dy<-pow(10,-10)||dz<-pow(10,-10)) return 0;

	ae=a[0]+(a[4]-a[0])*dz/GridSize0;
	dh=a[3]+(a[7]-a[3])*dz/GridSize0;
	bf=a[1]+(a[5]-a[1])*dz/GridSize0;
	cg=a[2]+(a[6]-a[2])*dz/GridSize0;

	ee=ae+(dh-ae)*dy/GridSize0;
	ff=bf+(cg-bf)*dy/GridSize0;

	*b=ee+(ff-ee)*dx/GridSize0;	//bn value
	*c=(ff-ee)/GridSize0;		//dbn2dx

	return 1;
}


void Trap_Runge_Kutta(double *Btx,double *Bty,double *Btz,double x0,double y0,double z0,double vx0,double vy0,double vz0,double m,double *x,double *y,double *z,double *vx,double *vy,double *vz,double dt)
{
//define fi(t,vx,x,vy,y,vz,z)
//this method has higher precision than Euler(Newtonian)

	double k11,k12,k13,k14;
	double k21,k22,k23,k24;
	double k31,k32,k33,k34;
	double k41,k42,k43,k44;
	double k51,k52,k53,k54;
	double k61,k62,k63,k64;

	double ax,ay,az;
	double xx,yy,zz;

	
		TrapAcceleration(Btx,Bty,Btz,x0,y0,z0,m,&ax,&ay,&az);
		k11=ax;//d(vx)/dt
		k21=vx0;//d(x)/dt
		k31=ay;//d(vy)/dt
		k41=vy0;//d(y)/dt
		k51=az;//d(vz)/dt
		k61=vz0;//d(z)/dt

		k22=vx0+dt/2*k11;	//medial vx
		xx=x0+dt/2*k21;		//medial x
		k42=vy0+dt/2*k31;	//medial vy
		yy=y0+dt/2*k41;		//medial y
		k62=vz0+dt/2*k51;	//medial vz
		zz=z0+dt/2*k61;		//medial z

		TrapAcceleration(Btx,Bty,Btz,xx,yy,zz,m,&ax,&ay,&az);
		k12=ax;
		k32=ay;
		k52=az;

		k23=vx0+dt/2*k12;
		xx=x0+dt/2*k22;	
		k43=vy0+dt/2*k32;	
		yy=y0+dt/2*k42;	
		k63=vz0+dt/2*k52;	
		zz=z0+dt/2*k62;	

		TrapAcceleration(Btx,Bty,Btz,xx,yy,zz,m,&ax,&ay,&az);
		k13=ax;
		k33=ay;
		k53=az;
				
		k24=vx0+dt*k13;	
		xx=x0+dt*k23;
		k44=vy0+dt*k33;	
		yy=y0+dt*k43;
		k64=vz0+dt*k53;	
		zz=z0+dt*k63;	

		TrapAcceleration(Btx,Bty,Btz,xx,yy,zz,m,&ax,&ay,&az);
		k14=ax;
		k34=ay;
		k54=az;

		*vx=vx0+dt/6*(k11+2*k12+2*k13+k14);
		*x=x0+dt/6*(k21+2*k22+2*k23+k24);	
		*vy=vy0+dt/6*(k31+2*k32+2*k33+k34);
		*y=y0+dt/6*(k41+2*k42+2*k43+k44);
		*vz=vz0+dt/6*(k51+2*k52+2*k53+k54);
		*z=z0+dt/6*(k61+2*k62+2*k63+k64);

}
