#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "allvars.h"
#include "proto.h"

/*! \file read_table.c
 *  \brief to get various kinds of tables calculated beforehand for different DEDMI models
 */
#ifdef HUBBLE_TABLE 
void me_init_hubble_table(void)
{
  FILE *fd;
  int loop;
  double ia,ih;
  double drift[HUBBLE_TABLE_LENGTH],gravkick[HUBBLE_TABLE_LENGTH],hydrokick[HUBBLE_TABLE_LENGTH];
  char fname[100];
  strcpy(fname,"hubble_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("counting Hubble Table length...\n");
	fflush(stdout);
  }  
  HUBBLE_TABLE_LENGTH=0;
  do
    {
      if(fscanf(fd, "%lg,%lg\n", &ia, &ih) == 2)
		{
			HUBBLE_TABLE_LENGTH++;
		}
      else
      {
		break;
		}
    }
  while(1);
  if(ThisTask==0)
  {
	printf("we find %d pairs of hubble table value.\n",HUBBLE_TABLE_LENGTH);
	fflush(stdout);
  }
  fclose(fd);
  double a[HUBBLE_TABLE_LENGTH],h[HUBBLE_TABLE_LENGTH];
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("reading `%s' ...\n", fname);
	fflush(stdout);
	printf("loading Hubble Table ...\n");
	fflush(stdout);
  } 
  for(loop=0;loop<HUBBLE_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&h[loop]);	
		//printf("input, a=%lf,h=%lf\n",a[loop],h[loop]); //debug use
		drift[loop] = 1.0/(All.Hubble*h[loop]*a[loop]*a[loop]*a[loop]);
		gravkick[loop] = 1.0/(All.Hubble*h[loop]*a[loop]*a[loop]);
		hydrokick[loop] = 1.0/(All.Hubble*h[loop]*pow(a[loop], 3 * GAMMA_MINUS1) * a[loop]);
	}
  MeHubbleAcc = gsl_interp_accel_alloc();
  MeHubbleSpline = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
  MeHubbleSplineDRIFT = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
  MeHubbleSplineGRAVKICK = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
  MeHubbleSplineHYDROKICK = gsl_spline_alloc(gsl_interp_cspline, HUBBLE_TABLE_LENGTH);
  gsl_spline_init(MeHubbleSpline, a, h, HUBBLE_TABLE_LENGTH);
  gsl_spline_init(MeHubbleSplineDRIFT, a, drift, HUBBLE_TABLE_LENGTH);
  gsl_spline_init(MeHubbleSplineGRAVKICK, a, gravkick, HUBBLE_TABLE_LENGTH);
  gsl_spline_init(MeHubbleSplineHYDROKICK, a, hydrokick, HUBBLE_TABLE_LENGTH);
  if(ThisTask==0)
  {
	printf("Hubble Table loading done.\n");
	fflush(stdout);
  }
  fclose(fd);
}
#endif
#ifdef DMMASS_TABLE
void me_init_dmmass_table(void)
{
  FILE *fd;
  int loop;
  char fname[100];
  double ia,im;
  strcpy(fname,"dmmass_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("counting DMMass Table length...\n");
	fflush(stdout);
  }  
  DMMASS_TABLE_LENGTH=0;
  do
    {
      if(fscanf(fd, "%lg,%lg\n", &ia, &im) == 2)
		{
			DMMASS_TABLE_LENGTH++;
		}
      else
      {
		break;
		}
    }
  while(1);
  if(ThisTask==0)
  {
	printf("we find %d pairs of dmmass table value.\n",DMMASS_TABLE_LENGTH);
	fclose(fd);
  }
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
    printf("loading DMMass Table ...\n");
    fflush(stdout);
  }  
  double a[DMMASS_TABLE_LENGTH],m[DMMASS_TABLE_LENGTH];
  for(loop=0;loop<DMMASS_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&m[loop]);
	}
  MeDMMassAcc = gsl_interp_accel_alloc();
  MeDMMassSpline = gsl_spline_alloc(gsl_interp_cspline, DMMASS_TABLE_LENGTH);
  gsl_spline_init(MeDMMassSpline, a, m, DMMASS_TABLE_LENGTH);
  if(ThisTask==0)
  {
    printf("DMMass Table loading done.\n");
    fflush(stdout);
  }
  fclose(fd);
}
#endif
#ifdef DRAG_TABLE
void me_init_drag_table(void)
{
  FILE *fd;
  int loop;
  char fname[100];
  double ia,id;
  strcpy(fname,"drag_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("counting Drag Table length...\n");
	fflush(stdout);
  }  
  DRAG_TABLE_LENGTH=0;
  do
    {
      if(fscanf(fd, "%lg,%lg\n", &ia, &id) == 2)
		{
			DRAG_TABLE_LENGTH++;
		}
      else
      {
		break;
		}
    }
  while(1);
  if(ThisTask==0)
  {
	printf("we find %d pairs of drag table value.\n",DRAG_TABLE_LENGTH);
	fclose(fd);
  }
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
    printf("loading Drag Table ...\n");
    fflush(stdout);
  }  
  double a[DRAG_TABLE_LENGTH],d[DRAG_TABLE_LENGTH];  
  for(loop=0;loop<DRAG_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&d[loop]);
	}
  All.MeDRAGAcc = gsl_interp_accel_alloc();
  All.MeDRAGSpline = gsl_spline_alloc(gsl_interp_cspline, DRAG_TABLE_LENGTH);
  gsl_spline_init(All.MeDRAGSpline, a, d, DRAG_TABLE_LENGTH);
  if(ThisTask==0)
  {
    printf("Drag Table loading done.\n");
    fflush(stdout);
  }
  fclose(fd);
}
#endif
#ifdef DE_TABLE
void me_init_de_table(void)
{
  printf("This function is not available in the public version yet.\n");
  printf("Please contact the author Jiajun Zhang for further support.\n");
}
#endif
#ifdef VARG_TABLE
void me_init_varg_table(void)
{
  FILE *fd;
  int loop;
  char fname[100];
  double ia,ig;
  strcpy(fname,"varg_table.txt");
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("counting Varg Table length...\n");
	fflush(stdout);
  }  
  VARG_TABLE_LENGTH=0;
  do
    {
      if(fscanf(fd, "%lg,%lg\n", &ia, &ig) == 2)
		{
			VARG_TABLE_LENGTH++;
		}
      else
      {
		break;
		}
    }
  while(1);
  if(ThisTask==0)
  {
	printf("we find %d pairs of varg table value.\n",VARG_TABLE_LENGTH);
	fclose(fd);
  }
  if(!(fd = fopen(fname, "r")))
	{
	  printf("can't open file `%s`\n", fname);
	  exit(0);
	}
  if(ThisTask==0)
  {
	printf("reading `%s' ...\n", fname);
    fflush(stdout);
    printf("loading Varg Table ...\n");
    fflush(stdout);
  }  
  double a[VARG_TABLE_LENGTH],g[VARG_TABLE_LENGTH];  
  for(loop=0;loop<VARG_TABLE_LENGTH;loop++)
    {
		fscanf(fd,"%lf,%lf\n",&a[loop],&g[loop]);
	}
  All.MeVARGAcc = gsl_interp_accel_alloc();
  All.MeVARGSpline = gsl_spline_alloc(gsl_interp_cspline, VARG_TABLE_LENGTH);
  gsl_spline_init(All.MeVARGSpline, a, g, VARG_TABLE_LENGTH);
  if(ThisTask==0)
  {
    printf("Varg Table loading done.\n");
    fflush(stdout);
  }
  fclose(fd);
}
#endif
