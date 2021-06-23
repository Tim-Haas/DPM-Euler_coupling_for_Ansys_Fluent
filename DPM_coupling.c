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
-x_momentum_source in source terms (Cell zone conditions->zone->source terms->x_momentum)
-y_momentum_source in source terms (Cell zone conditions->zone->source terms->y_momentum)
-z_momentum_source in source terms (Cell zone conditions->zone->source terms->z_momentum)


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
#define FLAG_INCLUDE_REMAINDER 0 /*in case the bubble size is similar to the grid spacing, roundoff errors of the weightings might become a problem. In case FLAG_INCLUDE_REMAINDER=1, the remainder is set as the weighting factor of the center cell so that the sum of weights is always 1. However, this might cause stability issues*/
#define COUPLING_SCHEME 2 /*0=forward coupling only, 1=backward coupling only (not recommended), 2=two-way coupling*/


/***********************************************************************************
	Auxillary functions
***********************************************************************************/

/*Deen_Polynom computes the share of momentum exchange of each affected cell*/
real Deen_Polynom(real xi,real xb,real nx)
	{
	/*Deen, N.G., M. van Sint Annaland, and J. Kuipers, Multi-scale modeling of dispersed
	  gas–liquid two-phase flow. Chemical engineering science, 2004. 59(8-9): p. 1853-1861*/
	/*xi is the Eulerian coordiante, xb is the Lagrangian coordinate, nx is the 
	  influencing sphere parameter*/
	real deenshare=0.0;
	deenshare=15./16.*((pow((xi-xb),4.)/pow(nx,5.))-(2.*pow((xi-xb),2.)/pow(nx,3.))+1./nx);
	return deenshare;
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
	real xc[ND_ND]; 
	C_CENTROID(xc,check,tc);
	int j=0;
	for (j=0;j<3;j++)
		{	
		if (fabs(xc[j]-P_POS(p)[j])>distance)
			{
			return 1;/*not influenced*/
			} 
		}
	return 0; /*influenced*/
	}


/**********************************************************************************
	Drag force auxillary functions
***********************************************************************************/
real Bozzano_drag (real Eotvos_no, real Reynolds_no, real liquid_viscosity, real liquid_density)
	{
	/*G. Bozzano and M. Dente, Computers & Chemical Engineering 2001, vol. 25, pp. 571-576*/
	real Mo=GRAVITY_1*liquid_viscosity*liquid_viscosity*liquid_viscosity*liquid_viscosity/(liquid_density*SURFACE_TENSION*SURFACE_TENSION*SURFACE_TENSION);
	real f_1=48./Reynolds_no*((1.+12.*pow(Mo, 0.333))/(1.+36.*pow(Mo,0.333)));
	real f_2=0.9*pow(Eotvos_no,1.5)/(1.4*(1.+30.*pow(Mo,1./6.))+pow(Eotvos_no,1.5));
	real f=f_1+f_2;
	real deform=(10.*(1.+1.3*pow(Mo,1./6.))+3.1*Eotvos_no)/(10.*(1.+1.3*pow(Mo,1./6.))+Eotvos_no);
	real cd=f*deform;
	return cd;
	}

	
real Dijkhuizen_drag (real Eotvos_no, real Reynolds_no)
	{
	/*W Dijkhuizen, I Roghair, M Van Sint Annaland and JAM Kuipers, Chemical Engineering Science 2010, vol. 65, pp. 1415-1426*/
	real mei=16./Reynolds_no*(1.+(2./(1.+16./Reynolds_no+3.315/sqrt(Reynolds_no))));
	real dijk=(4.*Eotvos_no)/(Eotvos_no+9.5);
	real cd=sqrt(mei*mei+dijk*dijk);
	return cd;
	}

	
/**********************************************************************************
	Lift force auxillary functions
***********************************************************************************/
real Tomiyama_lift(real Eotvos_no, real Reynolds_no, real particle_diam, real liquid_density, real particle_density)
	{
	/*Akio Tomiyama, Hidesada Tamai, Iztok Zun and Shigeo Hosokawa, Chemical Engineering Science 2002, vol. 57, pp. 1849-1858*/
	real cl=0.0;
	real dh=particle_diam*pow((1.+0.163*pow(Eotvos_no,0.757)),(1./3.));
	real eod=GRAVITY_1*(liquid_density-particle_density)*dh*dh/SURFACE_TENSION;
	real feod=0.00105*eod*eod*eod-0.0159*eod*eod-0.0204*eod+0.474;
	if(eod<=4)
		{
		cl=MIN(0.288*tanh(0.121*Reynolds_no),feod);
		}
	else if(eod>4. && eod<=10.)
		{
		cl=feod;
		}
	else if(eod>10.)
		{
		cl=-0.27;
		}
	return cl;	
	}

	
real Ziegenhein_lift(real Eotvos_no, real Reynolds_no, real particle_diam, real liquid_density, real particle_density)
	{
	/*T Ziegenhein, A Tomiyama and D Lucas, International Journal of Multiphase Flow 2018, vol. 108, pp. 11-24*/
	real cl=0.0;
	real dh=particle_diam*pow((1.+0.7*pow(Eotvos_no,0.7)),(1./3.));
	real eod=GRAVITY_1*(liquid_density-particle_density)*dh*dh/SURFACE_TENSION;
	real feod=0.5-0.1*eod+0.002*eod*eod;
	if(eod<=1.2)
		{
		cl=((Reynolds_no+16.)/(2.*(Reynolds_no+29.)));
		}
	else
		{
		cl=feod;
		}
	return cl;
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
		int no_marked_cells=0; 		/*number of marked cells*/
		int prev_no_marked_cells=0;	/*number of marked cells in previous cascade loop */
		int cascade_stage;			/*starting point for next cascade stage*/
		cell_t *cell_list=NULL;
		real *weightingfactor=NULL;
		cell_t c0, c1;
		Thread *t0, *t1, *tf;
		face_t f;
		real xc[ND_ND];
		real deenvalue[ND_ND]={0.0};
		real influencing_sphere=P_DIAM(p)*INFLUENCING_RADIUS_FACTOR/2.0;
		real sum=0.0;			/*sum of all weighting factors (should add up to 1.0)*/
		real remainder=0.0;		/*difference between weighting factors and 1.0, used to avoid missing data for coupling*/
		real mean_velocity[ND_ND]={0.0};
		real vorticity[ND_ND]={0.0};
		real relative_velocity[ND_ND]={0.0};
		real vort_cross_product[ND_ND]={0.0};
		real total_drag=0.0;
		real liftforce[ND_ND]={0.0};
		real dragforce[ND_ND]={0.0};
		real tot_relative_velocity=0.0;
		real cd=0.0;	/*Drag coefficient*/
		real cl=0.0;	/*Lift coefficient*/
		real Reynolds_no=0.0;
		real liquid_density=C_R(c,t);
		real liquid_viscosity=C_MU_L(c,t);
		real particle_diam=P_DIAM(p);
		real particle_density=P_RHO(p);
		real Eotvos_no=GRAVITY_1*(liquid_density-particle_density)*particle_diam*particle_diam/SURFACE_TENSION;
		real partmass=P_DIAM(p)*P_DIAM(p)*P_DIAM(p)/6.0*3.14*P_RHO(p);
		
		/* Test if influencing_sphere is smaller than cellsize, then whole impulse is stored in containing cell*/
		if (pow(C_VOLUME(c,t),1.0/3.0)>=influencing_sphere)
			{
			cell_list = (cell_t*) calloc(1,sizeof(cell_t));
			cell_list[no_marked_cells]=c;
			no_marked_cells++;
			prev_no_marked_cells=no_marked_cells;
			weightingfactor = (real*) calloc(1,sizeof(real));
			weightingfactor[0]=1.0;
			Message("single cell");
			}
			/* if influencing sphere is larger than cellsize start searching for cells*/
		else 
			{
			sum=0.0;
			/*initialize cell_list with cell the particle center is located*/
			cell_list = (cell_t*) calloc(1,sizeof(cell_t));
			cell_list[no_marked_cells]=c;
			no_marked_cells++;
			while (1) /*search for cells until no new cells are added*/
				{
				/*loop only over new cells because old cells have already chared their neighbors*/
				cascade_stage=prev_no_marked_cells;
				prev_no_marked_cells=no_marked_cells;
				for (arrayloop=cascade_stage; arrayloop<prev_no_marked_cells; arrayloop++)
					{
					c_face_loop(cell_list[arrayloop], t, n)         /* loops over all faces of a cell */
						{
						f = C_FACE(cell_list[arrayloop],t,n);
						tf = C_FACE_THREAD(cell_list[arrayloop],t,n);
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
								if (check_exist(c0,cell_list,no_marked_cells)==1)
									{ /*check if c0 is stored*/
									if (check_distance(c0,t,p,influencing_sphere)==0)
										{
										cell_list = (cell_t*) realloc(cell_list,(no_marked_cells+1)*sizeof(cell_t));
										if (cell_list != NULL)
											{
											cell_list[no_marked_cells]=c0;
											no_marked_cells++;
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
								if (check_exist(c1,cell_list,no_marked_cells)==1)
									{ /*repeat for c1*/
									if (check_distance(c1,t,p,influencing_sphere)==0)
										{ 
										cell_list = (cell_t*) realloc(cell_list,(no_marked_cells+1)*sizeof(cell_t));
										if (cell_list != NULL)
											{
											cell_list[no_marked_cells]=c1;
											no_marked_cells++;
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
				if (prev_no_marked_cells==no_marked_cells)
					{
					break;
					}
				}/*end of while endless loop*/	
		weightingfactor = (real*) calloc(no_marked_cells,sizeof(real));
		for (arrayloop=0; arrayloop<prev_no_marked_cells; arrayloop++)
			{
			if (C_VOF(cell_list[arrayloop], pt[PHASE_TO_DELETE])>= 0.1)
				{
				MARK_TP(p, P_FL_REMOVED);
				Message("Was deleted\n");
				return 0.;
				}
			C_CENTROID(xc,cell_list[arrayloop],t);
			for (dir=0;dir<3;dir++)
				{
				deenvalue[dir]=Deen_Polynom(xc[dir],P_POS(p)[dir],influencing_sphere);	
				}
			weightingfactor[arrayloop]=deenvalue[0]*deenvalue[1]*deenvalue[2]*C_VOLUME(cell_list[arrayloop],t);
			sum+=deenvalue[0]*deenvalue[1]*deenvalue[2]*C_VOLUME(cell_list[arrayloop],t);
			}
		} /*end of if-statement (smaller than influencing sphere)*/

	if (FLAG_INCLUDE_REMAINDER==1)
		{
		remainder=MAX((1.0-sum),0.0);
		if (remainder==0.0)
			{
			Message("Remainder zero\n");
			}
		weightingfactor[0]=remainder;}
		if (fabs(C_W(c,t))>10.0||fabs(P_VEL(p)[2])>10.0||fabs(C_W(c,t))<0.000000001||fabs(P_VEL(p)[2])<0.000000001)
			{
			for (dir=0;dir<ND_ND;dir++)
				{
				Message("Problem: CW:%lf PW:%lf\n",C_W(c,t),P_VEL(p)[2]);
				mean_velocity[dir]=0.5;
				vorticity[dir]=0.00001;
				P_VEL(p)[dir]=0.5;
				C_W(c,t)=0.5;
				}
			}
		else
			{
			if (cell_list != NULL && weightingfactor!=NULL)
				{
				for (arrayloop=0; arrayloop<prev_no_marked_cells; arrayloop++)
					{	
					mean_velocity[0]+=C_U(cell_list[arrayloop],t)*weightingfactor[arrayloop];
					mean_velocity[1]+=C_V(cell_list[arrayloop],t)*weightingfactor[arrayloop];
					mean_velocity[2]+=C_W(cell_list[arrayloop],t)*weightingfactor[arrayloop];	
					vorticity[0]+=(C_DWDY(cell_list[arrayloop],t)-C_DVDZ(cell_list[arrayloop],t))*weightingfactor[arrayloop];
					vorticity[1]+=(C_DUDZ(cell_list[arrayloop],t)-C_DWDX(cell_list[arrayloop],t))*weightingfactor[arrayloop];
					vorticity[2]+=(C_DVDX(cell_list[arrayloop],t)-C_DUDY(cell_list[arrayloop],t))*weightingfactor[arrayloop];
					}
				if (COUPLING_SCHEME==1)
					{
					mean_velocity[0]=C_U(c,t);
					mean_velocity[1]=C_V(c,t);
					mean_velocity[2]=C_W(c,t);
					vorticity[0]=(C_DWDY(c,t)-C_DVDZ(c,t));
					vorticity[1]=(C_DUDZ(c,t)-C_DWDX(c,t));
					vorticity[2]=(C_DVDX(c,t)-C_DUDY(c,t));
					}
				}
			} /*end of if-statement*/
	for (dir=0;dir<ND_ND;dir++)
		{
		relative_velocity[dir]=mean_velocity[dir]-P_VEL(p)[dir];
		}
	tot_relative_velocity=MAX(0.00001,sqrt(relative_velocity[0]*relative_velocity[0]+relative_velocity[1]*relative_velocity[1]+relative_velocity[2]*relative_velocity[2]));	 
	Reynolds_no=particle_diam*liquid_density*tot_relative_velocity/liquid_viscosity;	 
	if (DRAG_FLAG==1)
		{
		cd=MAX(MIN(16./Reynolds_no*(1.+0.15*pow(Reynolds_no,0.687)),48./Reynolds_no),8./3.*Eotvos_no/(Eotvos_no+4.)); /*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/
		}
	 else if (DRAG_FLAG==2)
		{
		 cd=MAX(MIN(24./Reynolds_no*(1.+0.15*pow(Reynolds_no,0.687)),72./Reynolds_no),8./3.*Eotvos_no/(Eotvos_no+4.));
		 /*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/
		}
	 else if (DRAG_FLAG==3)
		{
		cd=MAX(24./Reynolds_no*(1.+0.15*pow(Reynolds_no,0.687)),(8./3.)*Eotvos_no/(Eotvos_no+4.));
		/*Akio Tomiyama, Isao Kataoka, Iztok Zun and Tadashi Sakuguchi,
		JSME International Journal Series B Fluids and Thermal Engineering 1998, vol. 41, pp. 472-479*/ 
		}
	 else if (DRAG_FLAG==4)
		{
		cd=Bozzano_drag(Eotvos_no,Reynolds_no,liquid_viscosity,liquid_density);
		}
	else if (DRAG_FLAG==5)
		{
		cd=Dijkhuizen_drag (Eotvos_no,Reynolds_no);
		}
	else 
		{
		Message ("Error: No Drag defined");
		}
	if (LIFT_FLAG==1)
		{
		cl=Tomiyama_lift(Eotvos_no,Reynolds_no,particle_diam,liquid_density,particle_density);
		}
	else if(LIFT_FLAG==2) 
		{
		cl=Ziegenhein_lift(Eotvos_no,Reynolds_no,particle_diam,liquid_density,particle_density);
		}
	else if (LIFT_FLAG==3)
		{
		cl=0.5;
		}
	else
		{
		cl=0.0;
		}
	vort_cross_product[0]=relative_velocity[2]*vorticity[1]-relative_velocity[1]*vorticity[2];
	vort_cross_product[1]=relative_velocity[0]*vorticity[2]-relative_velocity[2]*vorticity[0];
	vort_cross_product[2]=relative_velocity[1]*vorticity[0]-relative_velocity[0]*vorticity[1];
	total_drag=0.5*cd*tot_relative_velocity*tot_relative_velocity*liquid_density*3.14*P_DIAM(p)*P_DIAM(p)/4.;

	for (dir=0;dir<ND_ND;dir++)
		{
		dragforce[dir]=total_drag*relative_velocity[dir]/tot_relative_velocity;
		liftforce[dir]=-1.0*cl*3.14*P_DIAM(p)*P_DIAM(p)*P_DIAM(p)/6.0*liquid_density*vort_cross_product[dir];
		P_USER_REAL(p,dir)=(dragforce[dir]+liftforce[dir])/partmass;
		}
	if (cell_list != NULL && weightingfactor!=NULL)
		{					
		for (arrayloop=0; arrayloop<prev_no_marked_cells; arrayloop++)
			{
			for (dir=0;dir<ND_ND;dir++) 
				{
				C_UDMI(cell_list[arrayloop],t,(dir))+=(-1.0)*(dragforce[dir]+liftforce[dir])*weightingfactor[arrayloop];
				if (COUPLING_SCHEME==0)
					{
					C_UDMI(c,t,(dir))+=((-1.0)*(dragforce[dir]+liftforce[dir]));
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
DEFINE_DPM_DRAG(zero_drag,Reynolds_no,p)
	{
	real drag_force;
	real cd=0.0;	
	drag_force=0.0;
    return (drag_force);    
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
	real x_momentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		x_momentum=MIN(MAX(C_UDMI(c,t,3)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return x_momentum;	
	#endif
	}


DEFINE_SOURCE(y_momentum_source,c,t,dS,eqn)
	{
	#if !RP_HOST	
	real y_momentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		y_momentum=MIN(MAX(C_UDMI(c,t,4)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return y_momentum;	
	#endif
	}


 DEFINE_SOURCE(z_momentum_source,c,t,dS,eqn)
	{
	#if !RP_HOST	
	real z_momentum=0.0;
	Thread *tt, **pt;
	Domain *domain=Get_Domain(MIXTURE_DOMAIN);
	tt=Lookup_Thread(domain,CELL_ZONE_ID_FLUID);
	pt = THREAD_SUB_THREADS(tt);
	if (C_VOF(c, pt[PHASE_TO_DELETE])< 0.25)
		{
		z_momentum=MIN(MAX(C_UDMI(c,t,5)/C_VOLUME(c,t),-1.0*GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP)),GRAVITY_1*(C_R(c,t)-GAS_DENSITY_STP));	
		}
	dS[eqn]=0;
	return z_momentum;	
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