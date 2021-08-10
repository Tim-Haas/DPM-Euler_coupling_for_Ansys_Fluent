/****************************************************************************
Coupling scheme for the Lagrangian DPM model in Ansys Fluent
Use for particles larger or similar to the Eulerian grid spacing. More realistic results 
by a fair share of the particle momentum to all covered cells. Makes the particles 
independent of the Eulerian grid, increases stability and prohibits pseudo-turbulence.
Steps:
1) Identifies all cells that are in the influencing sphere of a particle (bubble) by 
   an efficient cell loop cascade algorithm. 
2) Computes weighting factors for all covered cells based on the "Deen Polynom"
3) Maps the fluid flow properties to the particle
4) Maps the particle momentum to the Eulerian cells

Implemented with different drag and lift coefficients for bubbles
Compile UDF and set as follows:
-MappingMain in the DPM Model panel (Model->Discrete Phase->UDF->Scalar Update)
-ForwardCoupling in the DPM Model panel (Model->Discrete Phase->UDF->Body force)
-reset_particle_list in Execute at end (User-Defined->Function-hooks)
-zero_drag in DPM injection panel (Model->Discrete Phase->Injections->Physical Models->Drag Law)
-reset_udms in Execute at end (User-Defined->Function-hooks)
-x_momentum_source in source terms (Cell zone conditions->zone->source terms->xMomentum)
-y_momentum_source in source terms (Cell zone conditions->zone->source terms->yMomentum)
-z_momentum_source in source terms (Cell zone conditions->zone->source terms->zMomentum)


Written by: Tim Haas <haas@iob.rwth-aachen.de> and Christian Schubert in 2020
Department for Industrial Furnaces and Heat Engineering
RWTH Aachen University
*********************************************************************************/

/********************************************************************************
Copyright 2019-2020 Tim Haas
License (MIT):
Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in the 
Software without restriction, including without limitation the rights to use, copy, 
modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the 
following conditions:
The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ***********************************************************************************/

#include "udf.h"
#include "math.h"

/***********************************************************************************
	Global variables
you can customize the scheme to your application by adjusting the fluid properties
and choose between different drag and lift models
***********************************************************************************/

#define MIXTURE_DOMAIN 1
#define STD_PRESSURE_1 101325.0  		/*ambient pressure above the bath in Pa*/
#define BATH_HEIGTH_1 0.64				/*Filling height in m*/
#define SURFACE_TENSION 0.76
#define GAS_DENSITY_STP 1.2
#define GRAVITY_1 9.81
#define CELL_ZONE_ID_FLUID 7    		/*Get CELL_ZONE_ID_FLUID from the cell zone conditions Panel in fluent*/
#define INFLUENCING_RADIUS_FACTOR 4.0 	/*Radius of the sphere that indicates all affected Euerian cells, defined as a multiple of the bubble radius*/
#define PHASE_TO_DELETE 1    			/*ID of the VOF phase at which DPM particles should be deleted (first phase=0)*/
#define DRAG_FLAG 1 
/**********************************************************************************
	Implemented drag correlations:
1=Tomiyama drag correlation for pure systems
2=Tomiyama drag correlation for slightly contamined systems
3=Tomiyama drag correlation for fully contamined systems
4=Bozzano drag correlation
5=Dijkhuizen drag correlation
***********************************************************************************/
#define LIFT_FLAG 1
/**********************************************************************************
	Implemented lift correlations:
0=No lift force
1=Tomiyama lift correlation
2=Ziegenheim lift correlation
3=Const lift coefficient 0.5
***********************************************************************************/
#define VERSION 1.0
#define FLAG_INCLUDE_REMAINDER 0 /*in case the bubble size is similar to the grid spacing, roundoff errors of the weightings might become a problem. In case FLAG_INCLUDE_REMAINDER=1, the weightRemainder is set as the weighting factor of the center cell so that the sumOfWeights of weights is always 1. However, this might cause stability issues*/
#define COUPLING_SCHEME 2 /*0=forward coupling only, 1=backward coupling only (not recommended), 2=two-way coupling*/


/***********************************************************************************
	Auxillary functions
***********************************************************************************/

/*Deen_Polynom computes the share of momentum exchange of each affected cell*/
real Deen_Polynom(real eulerCoord,real lagrangeCoord,real inflSphere)
	{
	/*Deen, N.G., M. van Sint Annaland, and J. Kuipers, Multi-scale modeling of dispersed
	  gas–liquid two-phase flow. Chemical engineering science, 2004. 59(8-9): p. 1853-1861*/
	/*eulerCoord is the Eulerian coordiante, lagrangeCoord is the Lagrangian coordinate, inflSphere is the 
	  influencing sphere parameter*/
	real deenShare=0.0;
	deenShare=15./16.*((pow((eulerCoord-lagrangeCoord),4.)/pow(inflSphere,5.))-(2.*pow((eulerCoord-lagrangeCoord),2.)/pow(inflSphere,3.))+1./inflSphere);
	return deenShare;
	}

	
/*check_exist checks if a cell is already marked as affected (part of the cell array)*/
int check_exist (cell_t check, cell_t *carray, int length)
	{
	/*check is the cell thread of a cell that might be append to the cell list *carray*/ 
	int j=0;
	for (j=0; j<length; j++)
		{
		if (check==carray[j])
			{
			return 0;
			}
		}
	return 1;
	}


/*check_distance checks if a cell is within the influencing sphere*/
int check_distance (cell_t check, Thread *tc, Tracked_Particle *p, real distance)
	{
	real cellCenter[ND_ND]; 
	C_CENTROID(cellCenter,check,tc);
	int j=0;
	for (j=0;j<3;j++)
		{	
		if (fabs(cellCenter[j]-P_POS(p)[j])>distance)
			{
			return 1;/*not influenced*/
			} 
		}
	return 0; /*influenced*/
	}


/**********************************************************************************
	Drag force auxillary functions
***********************************************************************************/
real Bozzano_drag (real EotvosNo, real ReynoldsNo, real liquidViscosity, real liquidDensity)
	{
	/*G. Bozzano and M. Dente, Computers & Chemical Engineering 2001, vol. 25, pp. 571-576*/
	real Morton_no=GRAVITY_1*liquidViscosity*liquidViscosity*liquidViscosity*liquidViscosity/(liquidDensity*SURFACE_TENSION*SURFACE_TENSION*SURFACE_TENSION);
	real fOne=48./ReynoldsNo*((1.+12.*pow(Morton_no, 0.333))/(1.+36.*pow(Morton_no,0.333)));
	real fTwo=0.9*pow(EotvosNo,1.5)/(1.4*(1.+30.*pow(Morton_no,1./6.))+pow(EotvosNo,1.5));
	real f=fOne+fTwo;
	real deform=(10.*(1.+1.3*pow(Morton_no,1./6.))+3.1*EotvosNo)/(10.*(1.+1.3*pow(Morton_no,1./6.))+EotvosNo);
	real dragCoef=f*deform;
	return dragCoef;
	}

	
real Dijkhuizen_drag (real EotvosNo, real ReynoldsNo)
	{
	/*W Dijkhuizen, I Roghair, M Van Sint Annaland and JAM Kuipers, Chemical Engineering Science 2010, vol. 65, pp. 1415-1426*/
	real mei=16./ReynoldsNo*(1.+(2./(1.+16./ReynoldsNo+3.315/sqrt(ReynoldsNo))));
	real dijk=(4.*EotvosNo)/(EotvosNo+9.5);
	real dragCoef=sqrt(mei*mei+dijk*dijk);
	return dragCoef;
	}

	
/**********************************************************************************
	Lift force auxillary functions
***********************************************************************************/
real Tomiyama_lift(real EotvosNo, real ReynoldsNo, real particleDiam, real liquidDensity, real particleDensity)
	{
	/*Akio Tomiyama, Hidesada Tamai, Iztok Zun and Shigeo Hosokawa, Chemical Engineering Science 2002, vol. 57, pp. 1849-1858*/
	real liftCoef=0.0;
	real majAxis=particleDiam*pow((1.+0.163*pow(EotvosNo,0.757)),(1./3.));
	real defEotvos=GRAVITY_1*(liquidDensity-particleDensity)*majAxis*majAxis/SURFACE_TENSION;
	real feod=0.00105*defEotvos*defEotvos*defEotvos-0.0159*defEotvos*defEotvos-0.0204*defEotvos+0.474;
	if(defEotvos<=4)
		{
		liftCoef=MIN(0.288*tanh(0.121*ReynoldsNo),feod);
		}
	else if(defEotvos>4. && defEotvos<=10.)
		{
		liftCoef=feod;
		}
	else if(defEotvos>10.)
		{
		liftCoef=-0.27;
		}
	return liftCoef;	
	}

	
real Ziegenhein_lift(real EotvosNo, real ReynoldsNo, real particleDiam, real liquidDensity, real particleDensity)
	{
	/*T Ziegenhein, A Tomiyama and D Lucas, International Journal of Multiphase Flow 2018, vol. 108, pp. 11-24*/
	real liftCoef=0.0;
	real majAxis=particleDiam*pow((1.+0.7*pow(EotvosNo,0.7)),(1./3.));
	real defEotvos=GRAVITY_1*(liquidDensity-particleDensity)*majAxis*majAxis/SURFACE_TENSION;
	real feod=0.5-0.1*defEotvos+0.002*defEotvos*defEotvos;
	if(defEotvos<=1.2)
		{
		liftCoef=((ReynoldsNo+16.)/(2.*(ReynoldsNo+29.)));
		}
	else
		{
		liftCoef=feod;
		}
	return liftCoef;
	}


/***********************************************************************************
	Main function
	Finds all cells that are covered by the influencing sphere of a bubble
	Computes the force acting on the bubble by using the weighted flow variables of
	Eulerian cells and store them in a Particle Variable (P_USER_REAL)
	Compute the exchange momentum for each Eulerian cell and store them in a 
	User Defined Memory Variable (UDMI)
***********************************************************************************/
DEFINE_DPM_SCALAR_UPDATE(MappingMain, c, t, init, p)
	{
	#if !RP_HOST
	
	/*Definition of a pointer to Cell_Thread to determine whether the bubble has reached the VOF phase boundary (e.g. free surface)*/
	Thread *tt, **pt;
	int64_t *realloc_particlelist = NULL;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);

	/*Bubble grow*/	
	P_DIAM(p) = P_INIT_DIAM(p)*pow((C_R(c,t)*GRAVITY_1*(BATH_HEIGTH_1-P_INIT_POS(p)[2])+STD_PRESSURE_1)/(C_R(c,t)*GRAVITY_1*(BATH_HEIGTH_1-P_POS(p)[2])+STD_PRESSURE_1), 1/3.);
	P_RHO(p)=P_INIT_RHO(p)*(STD_PRESSURE_1+C_R(c,t)*GRAVITY_1*(BATH_HEIGTH_1-P_POS(p)[2]))/(STD_PRESSURE_1+C_R(c,t)*GRAVITY_1*(BATH_HEIGTH_1-P_INIT_POS(p)[2]));
	
	/*Delete criterion at phase boundary*/
	if (C_VOF(c, pt[PHASE_TO_DELETE])>= 0.1)
		{
		MARK_TP(p, P_FL_REMOVED);
		return 0.;
		}
		
	if (init == 1) /*DEFINE_DPM_SCALAR_UPDATE is called twice by Fluent. The init variable indicates whether it is called at the beginning or end of a time step*/
		{	
		int dir, n, arrayloop; 		/*loop variables*/
		int noMarkedCells=0; 		/*number of marked cells*/
		int prevNoMarkedCells=0;	/*number of marked cells in previous cascade loop */
		int cascadeStage;			/*starting point for next cascade stage*/
		cell_t *cellList=NULL;
		real *weightingFactor=NULL;
		cell_t c0, c1;
		Thread *t0, *t1, *tf;
		face_t f;
		real xc[ND_ND];
		real deenValue[ND_ND]={0.0};
		real influencingSphere=P_DIAM(p)*INFLUENCING_RADIUS_FACTOR/2.0;
		real sumOfWeights=0.0;			/*sumOfWeights of all weighting factors (should add up to 1.0)*/
		real weightRemainder=0.0;		/*difference between weighting factors and 1.0, used to avoid missing data for coupling*/
		real meanVelocity[ND_ND]={0.0};
		real vorticity[ND_ND]={0.0};
		real relativeVelocity[ND_ND]={0.0};
		real vortCrossProduct[ND_ND]={0.0};
		real totalDrag=0.0;
		real liftForce[ND_ND]={0.0};
		real dragForce[ND_ND]={0.0};
		real totRelativeVelocity=0.0;
		real dragCoef=0.0;	/*Drag coefficient*/
		real liftCoef=0.0;	/*Lift coefficient*/
		real ReynoldsNo=0.0;
		real liquidDensity=C_R(c,t);
		real liquidViscosity=C_MU_L(c,t);
		real particleDiam=P_DIAM(p);
		real particleDensity=P_RHO(p);
		real EotvosNo=GRAVITY_1*(liquidDensity-particleDensity)*particleDiam*particleDiam/SURFACE_TENSION;
		real particleMass=P_DIAM(p)*P_DIAM(p)*P_DIAM(p)/6.0*3.14*P_RHO(p);
		
		/* Test if influencingSphere is smaller than cellsize, then whole impulse is stored in containing cell*/
		if (pow(C_VOLUME(c,t),1.0/3.0)>=influencingSphere)
			{
			cellList = (cell_t*) calloc(1,sizeof(cell_t));
			cellList[noMarkedCells]=c;
			noMarkedCells++;
			prevNoMarkedCells=noMarkedCells;
			weightingFactor = (real*) calloc(1,sizeof(real));
			weightingFactor[0]=1.0;
			Message("single cell");
			}
			/* if influencing sphere is larger than cellsize start searching for cells*/
		else 
			{
			sumOfWeights=0.0;
			/*initialize cellList with cell the particle center is located*/
			cellList = (cell_t*) calloc(1,sizeof(cell_t));
			cellList[noMarkedCells]=c;
			noMarkedCells++;
			while (1) /*search for cells until no new cells are added*/
				{
				/*loop only over new cells because old cells have already chared their neighbors*/
				cascadeStage=prevNoMarkedCells;
				prevNoMarkedCells=noMarkedCells;
				for (arrayloop=cascadeStage; arrayloop<prevNoMarkedCells; arrayloop++)
					{
					c_face_loop(cellList[arrayloop], t, n)         /* loops over all faces of a cell */
						{
						f = C_FACE(cellList[arrayloop],t,n);
						tf = C_FACE_THREAD(cellList[arrayloop],t,n);
						if(BOUNDARY_FACE_THREAD_P(tf))
							{
							continue;
							}
						else
							{
							c0=F_C0(f,tf);
							t0=THREAD_T0(tf);
							c1=F_C1(f,tf);
							t1=THREAD_T1(tf);
							if (t0==t)
								{
								/*loop over stored cells if cell id is already stored*/
								if (check_exist(c0,cellList,noMarkedCells)==1)
									{ /*check if c0 is stored*/
									if (check_distance(c0,t,p,influencingSphere)==0)
										{
										cellList = (cell_t*) realloc(cellList,(noMarkedCells+1)*sizeof(cell_t));
										if (cellList != NULL)
											{
											cellList[noMarkedCells]=c0;
											noMarkedCells++;
											}
										else
											{
											Message("Memory realloc error!\n");
											}
										}								
									}
								}
							if (t1==t)
								{
								if (check_exist(c1,cellList,noMarkedCells)==1)
									{ /*repeat for c1*/
									if (check_distance(c1,t,p,influencingSphere)==0)
										{ 
										cellList = (cell_t*) realloc(cellList,(noMarkedCells+1)*sizeof(cell_t));
										if (cellList != NULL)
											{
											cellList[noMarkedCells]=c1;
											noMarkedCells++;
											}
										else
											{
											Message("Memory realloc error!\n");
											}					
										}			
									}
								}
							}
						}/*end of face loop*/
					} /*end of loop over cells*/
				if (prevNoMarkedCells==noMarkedCells)
					{
					break;
					}
				}/*end of while endless loop*/	
		weightingFactor = (real*) calloc(noMarkedCells,sizeof(real));
		for (arrayloop=0; arrayloop<prevNoMarkedCells; arrayloop++)
			{
			if (C_VOF(cellList[arrayloop], pt[PHASE_TO_DELETE])>= 0.1)
				{
				MARK_TP(p, P_FL_REMOVED);
				Message("Was deleted\n");
				return 0.;
				}
			C_CENTROID(xc,cellList[arrayloop],t);
			for (dir=0;dir<3;dir++)
				{
				deenValue[dir]=Deen_Polynom(xc[dir],P_POS(p)[dir],influencingSphere);	
				}
			weightingFactor[arrayloop]=deenValue[0]*deenValue[1]*deenValue[2]*C_VOLUME(cellList[arrayloop],t);
			sumOfWeights+=deenValue[0]*deenValue[1]*deenValue[2]*C_VOLUME(cellList[arrayloop],t);
			}
		} /*end of if-statement (smaller than influencing sphere)*/

	if (FLAG_INCLUDE_REMAINDER==1)
		{
		weightRemainder=MAX((1.0-sumOfWeights),0.0);
		if (weightRemainder==0.0)
			{
			Message("weightRemainder zero\n");
			}
		weightingFactor[0]=weightRemainder;}
		if (fabs(C_W(c,t))>10.0||fabs(P_VEL(p)[2])>10.0||fabs(C_W(c,t))<0.000000001||fabs(P_VEL(p)[2])<0.000000001)
			{
			for (dir=0;dir<ND_ND;dir++)
				{
				Message("Problem: CW:%lf PW:%lf\n",C_W(c,t),P_VEL(p)[2]);
				meanVelocity[dir]=0.5;
				vorticity[dir]=0.00001;
				P_VEL(p)[dir]=0.5;
				C_W(c,t)=0.5;
				}
			}
		else
			{
			if (cellList != NULL && weightingFactor!=NULL)
				{
				for (arrayloop=0; arrayloop<prevNoMarkedCells; arrayloop++)
					{	
					meanVelocity[0]+=C_U(cellList[arrayloop],t)*weightingFactor[arrayloop];
					meanVelocity[1]+=C_V(cellList[arrayloop],t)*weightingFactor[arrayloop];
					meanVelocity[2]+=C_W(cellList[arrayloop],t)*weightingFactor[arrayloop];	
					vorticity[0]+=(C_DWDY(cellList[arrayloop],t)-C_DVDZ(cellList[arrayloop],t))*weightingFactor[arrayloop];
					vorticity[1]+=(C_DUDZ(cellList[arrayloop],t)-C_DWDX(cellList[arrayloop],t))*weightingFactor[arrayloop];
					vorticity[2]+=(C_DVDX(cellList[arrayloop],t)-C_DUDY(cellList[arrayloop],t))*weightingFactor[arrayloop];
					}
				if (COUPLING_SCHEME==1)
					{
					meanVelocity[0]=C_U(c,t);
					meanVelocity[1]=C_V(c,t);
					meanVelocity[2]=C_W(c,t);
					vorticity[0]=(C_DWDY(c,t)-C_DVDZ(c,t));
					vorticity[1]=(C_DUDZ(c,t)-C_DWDX(c,t));
					vorticity[2]=(C_DVDX(c,t)-C_DUDY(c,t));
					}
				}
			} /*end of if-statement*/
	for (dir=0;dir<ND_ND;dir++)
		{
		relativeVelocity[dir]=meanVelocity[dir]-P_VEL(p)[dir];
		}
	totRelativeVelocity=MAX(0.00001,sqrt(relativeVelocity[0]*relativeVelocity[0]+relativeVelocity[1]*relativeVelocity[1]+relativeVelocity[2]*relativeVelocity[2]));	 
	ReynoldsNo=particleDiam*liquidDensity*totRelativeVelocity/liquidViscosity;	 
	if (DRAG_FLAG==1)
		{
		dragCoef=MAX(MIN(16./ReynoldsNo*(1.+0.15*pow(ReynoldsNo,0.687)),48./ReynoldsNo),8./3.*EotvosNo/(EotvosNo+4.)); /*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/
		}
	 else if (DRAG_FLAG==2)
		{
		 dragCoef=MAX(MIN(24./ReynoldsNo*(1.+0.15*pow(ReynoldsNo,0.687)),72./ReynoldsNo),8./3.*EotvosNo/(EotvosNo+4.));
		 /*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/
		}
	 else if (DRAG_FLAG==3)
		{
		dragCoef=MAX(24./ReynoldsNo*(1.+0.15*pow(ReynoldsNo,0.687)),(8./3.)*EotvosNo/(EotvosNo+4.));
		/*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/ 
		}
	 else if (DRAG_FLAG==4)
		{
		dragCoef=Bozzano_drag(EotvosNo,ReynoldsNo,liquidViscosity,liquidDensity);
		}
	else if (DRAG_FLAG==5)
		{
		dragCoef=Dijkhuizen_drag (EotvosNo,ReynoldsNo);
		}
	else 
		{
		Message ("Error: No Drag defined");
		}
	if (LIFT_FLAG==1)
		{
		liftCoef=Tomiyama_lift(EotvosNo,ReynoldsNo,particleDiam,liquidDensity,particleDensity);
		}
	else if(LIFT_FLAG==2) 
		{
		liftCoef=Ziegenhein_lift(EotvosNo,ReynoldsNo,particleDiam,liquidDensity,particleDensity);
		}
	else if (LIFT_FLAG==3)
		{
		liftCoef=0.5;
		}
	else
		{
		liftCoef=0.0;
		}
	vortCrossProduct[0]=relativeVelocity[2]*vorticity[1]-relativeVelocity[1]*vorticity[2];
	vortCrossProduct[1]=relativeVelocity[0]*vorticity[2]-relativeVelocity[2]*vorticity[0];
	vortCrossProduct[2]=relativeVelocity[1]*vorticity[0]-relativeVelocity[0]*vorticity[1];
	totalDrag=0.5*dragCoef*totRelativeVelocity*totRelativeVelocity*liquidDensity*3.14*P_DIAM(p)*P_DIAM(p)/4.;

	for (dir=0;dir<ND_ND;dir++)
		{
		dragForce[dir]=totalDrag*relativeVelocity[dir]/totRelativeVelocity;
		liftForce[dir]=-1.0*liftCoef*3.14*P_DIAM(p)*P_DIAM(p)*P_DIAM(p)/6.0*liquidDensity*vortCrossProduct[dir];
		P_USER_REAL(p,dir)=(dragForce[dir]+liftForce[dir])/particleMass;
		}
	if (cellList != NULL && weightingFactor!=NULL)
		{					
		for (arrayloop=0; arrayloop<prevNoMarkedCells; arrayloop++)
			{
			for (dir=0;dir<ND_ND;dir++) 
				{
				C_UDMI(cellList[arrayloop],t,(dir))+=(-1.0)*(dragForce[dir]+liftForce[dir])*weightingFactor[arrayloop];
				if (COUPLING_SCHEME==0)
					{
					C_UDMI(c,t,(dir))+=((-1.0)*(dragForce[dir]+liftForce[dir]));
					}
				}
			}	
		}
	}/*end of init=1 statement*/
	#endif
	}


/***********************************************************************************
	Introduce the computed body forces (drag and lift) to the DPM force balance
***********************************************************************************/
DEFINE_DPM_BODY_FORCE(ForwardCoupling,p,i)
	{
	#if !RP_HOST
    return (P_USER_REAL(p,i));
	# endif
	}
	
	
/***********************************************************************************
	Set the drag coefficient in the default model to zero to avoid considering the 
	drag twice
***********************************************************************************/
DEFINE_DPM_DRAG(zero_drag,ReynoldsNo,p)
	{
	real dragForce;
	real dragCoef=0.0;	
	dragForce=0.0;
    return (dragForce);    
	}

	
/***********************************************************************************
	Reset the User Defined Memories for the next time step
***********************************************************************************/

DEFINE_EXECUTE_AT_END(reset_udms)
	{
	#if !RP_HOST
	Domain *domain=Get_Domain(1);
	cell_t c;
	Thread *t;
	thread_loop_c (t,domain)
		{
		begin_c_loop (c,t)
			{
			C_UDMI(c,t,3)=C_UDMI(c,t,0);	
			C_UDMI(c,t,4)=C_UDMI(c,t,1);	
			C_UDMI(c,t,5)=C_UDMI(c,t,2);	
			C_UDMI(c,t,0)=0.0;
			C_UDMI(c,t,1)=0.0;
			C_UDMI(c,t,2)=0.0;
			}
       end_c_loop (c,t)   
	   }
	#endif
	}
 
 
 /***********************************************************************************
	Introduce the exchanged momentum in the Navier Stokes Equations of the 
	Eulerian cells
***********************************************************************************/
DEFINE_SOURCE(x_momentum_source,c,t,dS,eqn)
	{
	#if !RP_HOST	
	real xMomentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		xMomentum=MIN(MAX(C_UDMI(c,t,3)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return xMomentum;	
	#endif
	}


DEFINE_SOURCE(y_momentum_source,c,t,dS,eqn)
	{
	#if !RP_HOST	
	real yMomentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		yMomentum=MIN(MAX(C_UDMI(c,t,4)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return yMomentum;	
	#endif
	}


 DEFINE_SOURCE(z_momentum_source,c,t,dS,eqn)
	{
	#if !RP_HOST	
	real zMomentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		zMomentum=MIN(MAX(C_UDMI(c,t,5)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return zMomentum;	
	#endif
	}
	
	
/***********************************************************************************
	Initialize the User Defined Memories
***********************************************************************************/
DEFINE_EXECUTE_ON_LOADING(initialize, libname)
	{
	#if !RP_HOST
	Domain *domain=Get_Domain(1);
	cell_t c;
	Thread *t;
	thread_loop_c (t,domain)
		{
		begin_c_loop (c,t)
			{
			C_UDMI(c,t,0)=0.0;
			C_UDMI(c,t,1)=0.0;
			C_UDMI(c,t,2)=0.0;
			C_UDMI(c,t,3)=0.0;
			C_UDMI(c,t,4)=0.0;	
			C_UDMI(c,t,5)=0.0;
			}
       end_c_loop (c,t)   
		}
	Message("Initialize Memory v.%lf\n",VERSION);
	#endif
	}
