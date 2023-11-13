// This routine writes out a 2D chop from the 3D solution
// Requires a list of nodes to create a chop X-section from
#include <stdlib.h>
#include "stdio.h"
#include "math.h"
#define MAXLINE 1000
#define MAXCARD 3
#define NDPRELM 4
#define CARD_E4T 8012
#define CARD_ND  19314
#define CARD_NDS 19333
#define UNSET_INT -3 /* an integer that has not been set */

/* Prototype Functions */
int init_split(
  char *,
  char *
);

int init_read_idatum(
  char *,
  char *
);

double init_read_ddatum(
  char *,
  char *
);
 
void init_read_error(
  char *,
  char *
); 
/* End Prototype definitions */

/* STRUCTURES */

typedef struct{
  double x, y, z;
  } VECT;
  
  
void tl_3d_to_2d()
{

/* START MAIN */
int i,j,k,l,m,n;
int nnode;
int nelem3d;
int isvel, isdpl,iscon;
int no_con;
int number;
int number3;
int *closest;
int n_tot;
int count;
double x_a, y_a, z_a;
double xx, yy, zz;
double *distance;
double TS, dp;
double vx,vy,vz;
double co;
double limit;
double scl;
FILE *input;
FILE *disp;
FILE *vel;
FILE *con;
FILE *mesh;
FILE *out_dpl;
FILE *out_vel;
FILE *out_con;
FILE *twodm;
FILE *nfile;
char line[MAXLINE];
char data[MAXLINE];
char str[50];
char name[50];
char str2[50] = ".3dm";
char str3[5] = ".dat";
char str4[50] = ".2dm";
char str5[50] = "mesh2d";
char str6[5] = ".nds";
char dpl[50] = "_dpl";
char vele[5] = "_vel";
char conc[6] = "_con";

VECT *node;
VECT *to_anal;
nelem3d = 0;
n_tot = 0;
count = 0;
input = fopen("CHOP_3D_ADH.txt","r");

/* READ DATA FFROM INPUT FILE */
while (!feof(input))
 {
   //fscanf(input,"%s\n",&str);
   fscanf(input,"%s\n",str);
   fscanf(input,"%d %d %d %d\n",&isdpl, &isvel, &iscon, &no_con);
   fscanf(input,"%lf\n",&scl);


 }


 /* OPEN MESH .3dm file */
 sprintf(name,"%s%s",str,str2);
 printf("\nThe ADH 3DM file is %s\n",name);
 mesh=fopen(name,"r");
 sprintf(name,"%s%s",str,str6);
 printf("The Chop Node List file is %s\n",name);
 nfile = fopen(name,"r");



while (fgets (line, MAXLINE, mesh) != NULL)
    {
      /* splits the line and goes to the appropriate case */
      switch (init_split (line, data))
	{
	case CARD_E4T:
	  number3 = init_read_idatum (line, data);

	  if (number3 > nelem3d)
	    nelem3d = number3;
	  break;
	default:
	  break;
	}
    }
  rewind (mesh);
  
  node = (VECT *) malloc (sizeof (VECT)* nelem3d * 4);
  nnode = 0;
 while (fgets (line, MAXLINE, mesh) != NULL)
    {
      /* splits the line and goes to the appropriate case */
      switch (init_split (line, data))
	{
	case CARD_ND:
      number = init_read_idatum (line, data);
	  xx = init_read_ddatum (line, data);
	  yy = init_read_ddatum (line, data);
	  zz = init_read_ddatum (line, data);
	  node[number - 1].x = xx;
	  node[number - 1].y = yy;
	  node[number - 1].z = zz;
	if (number > nnode)
	   nnode = number;
	  break;
	default:
	  break;
	}
    }
  rewind (mesh);
  printf("Total number of nodes in file is %d\n",nnode);
 number3 = 0;
while (fgets (line, MAXLINE, nfile) != NULL)
    {
      /* splits the line and goes to the appropriate case */
      switch (init_split (line, data))
	{
	case CARD_NDS:
	  number3 += 1;
	  break;
	default:
	  break;
	}
    }
printf("Total Nodes in Node String file are %d\n",number3);
rewind(nfile);	


to_anal = (VECT *) malloc (sizeof (VECT)* number3);
n_tot = number3;
number3 = 0;
i= 0;
while (fgets (line, MAXLINE, nfile) != NULL)
    {
      /* splits the line and goes to the appropriate case */
      switch (init_split (line, data))
	{
	case CARD_NDS:
	  { 
	  number3 = init_read_idatum (line, data);
	  to_anal[i].x = node[number3 -1].x;
	  to_anal[i].y = node[number3 -1].y;
	  i += 1;
	  }
	  break;
	default:
	  break;
	}
    }
	
// for (i =0; i < n_tot; i ++)
//  printf("X and Y for nodes to be analyzed are %lf %lf\n",to_anal[i].x, to_anal[i].y);

	   
	      
   printf("\nWriting 2D Chop\n");
    sprintf(name,"%s_CHOP%s",str, str4);
   twodm = fopen(name,"w");  
   fprintf(twodm,"mesh2d\n");
   for (i = 0;i < nnode; i ++)
   for (j = 0; j < n_tot;  j++)
	   if (node[i].x == to_anal[j].x && node[i].y == to_anal[j].y)
	   {
	   count +=1;
	    fprintf(twodm,"ND %d %lf %lf %lf\n",count,node[i].x, node[i].z * scl, node[i].y);
	   }
	   fflush(twodm);
	   fclose(twodm);

   if (iscon == 1)
     for (k = 0; k < no_con; k ++) 
{
     sprintf(name,"%s_con%d_chop%s",str,k+1,str3);
     out_con = fopen(name,"w");
	 fprintf(out_con,"DATASET\n");
	 fprintf(out_con,"OBJTYPE %s\n",str5);
     fprintf(out_con,"BEGSCL\n");
     fprintf(out_con,"ND %d\n",count);
	 fprintf(out_con,"NAME icon%d\n", k+1);
     sprintf(name,"%s_con%d%s",str,k+1,str3);
     con = fopen(name,"r");	
	 for(i = 0; i < 6; i ++)
     fgets(line,MAXLINE,con);
{ 
     while(!feof(con))
  {
    fscanf(con,"%*s %*d %lf\n",&TS);
	fprintf(out_con,"TS 0 %lf\n",TS);
		    
	  for (i = 0; i < nnode; i ++)
	  {
	  fscanf(con,"%lf\n",&co);
	    for (j = 0; j < n_tot; j ++)
	   {
		 if (node[i].x == to_anal[j].x && node[i].y == to_anal[j].y)
	        fprintf(out_con,"%lf\n",co);
         fflush(out_con);			
	   }
	   
	   }
   }
   fprintf(out_con,"ENDDS");
   } 


fclose(out_con);
fclose(con);
	 }


if (isvel ==1)
{
     sprintf(name,"%s_vel_chop%s",str,str3);
     out_vel = fopen(name,"w");
	 fprintf(out_vel,"DATASET\n");
	 fprintf(out_vel,"OBJTYPE %s\n",str5);
     fprintf(out_vel,"BEGVEC\n");
     fprintf(out_vel,"ND %d\n",count);
	 fprintf(out_vel,"NAME Velocity\n");
     sprintf(name,"%s_vel%s",str,str3);
     vel = fopen(name,"r");
	 
	 for(i = 0; i < 7; i ++)
     fgets(line,MAXLINE,vel);
{ 
     while(!feof(vel))
  {
    fscanf(vel,"%*s %*d %lf\n",&TS);
	fprintf(out_vel,"TS 0 %lf\n",TS);
		    
	  for (i = 0; i < nnode; i ++)
	  {
	  fscanf(vel,"%lf %lf %lf\n",&vx, &vy, &vz);
	    for (j = 0; j < n_tot; j ++)
	   {
		 if (node[i].x == to_anal[j].x && node[i].y == to_anal[j].y)
	        fprintf(out_vel,"%lf %lf %lf\n",vx,vz,vy);
         fflush(out_vel);			
	   }
	   
	   }
   }
   fprintf(out_vel,"ENDDS");
   } 
	 
	 
}	
fclose(out_vel);
fclose(vel); 

}	 



 int init_split(
  char *line,			/* the input line */
  char *data			/* the data */
)
{
  int ipmc;			/* i + MAXCARD */
  int card;			/* the card value */
  int card1, card2, card3;	/* temporary holders for the three card values */

  /* reads the card */
  /* for each card entry we assign it a number based on the following:
     (' ':0), (a-z:1-26), (A-Z:1-26), (0-9:27-36)

     the integer card value is given by:
     1369*i1+37*i2+i3
   */
  card1 = line[0];
  card2 = line[1];
  card3 = line[2];
  if(card1 >= 'a' && card1 <= 'z')
    card1 = card1 - 'a' + 1;
  else if(card1 >= 'A' && card1 <= 'Z')
    card1 = card1 - 'A' + 1;
  else if(card1 >= '0' && card1 <= '9')
    card1 = card1 - '0' + 27;
  else
    card1 = 0;
  if(card2 >= 'a' && card2 <= 'z')
    card2 = card2 - 'a' + 1;
  else if(card2 >= 'A' && card2 <= 'Z')
    card2 = card2 - 'A' + 1;
  else if(card2 >= '0' && card2 <= '9')
    card2 = card2 - '0' + 27;
  else
    card2 = 0;
  if(card3 >= 'a' && card3 <= 'z')
    card3 = card3 - 'a' + 1;
  else if(card3 >= 'A' && card3 <= 'Z')
    card3 = card3 - 'A' + 1;
  else if(card3 >= '0' && card3 <= '9')
    card3 = card3 - '0' + 27;
  else
    card3 = 0;
  card = 1369 * card1 + 37 * card2 + card3;

  /* sets the data array */
  for(ipmc = MAXCARD; ipmc < MAXLINE; ipmc++)
    {
      if(line[ipmc] == '!')
	{
	  data[ipmc - MAXCARD] = '\0';
	  return (card);
	}
      else if(line[ipmc] == '\0')
	{
	  data[ipmc - MAXCARD] = '\0';
	  return (card);
	}
      else
	{
	  data[ipmc - MAXCARD] = line[ipmc];
	}
    }

  /* if it drops out the bottom, then the line is too long */
 /* init_read_error(line, "Bad input line in init_split - too long.");*/
  return (UNSET_INT);
}


double init_read_ddatum(
  char *line,			/* the original line */
  char *data			/* the data */
)
{
  int i, j;			/* loop counter */
  int icnt;			/* white space */
  char datum[MAXLINE];		/* holds the datum */
  double value;			/* holds the value */

  /* initializes the datum */
  for(i = 0; i < MAXLINE; i++)
    datum[i] = '\0';

  /* strips the initial white space */
  icnt = 0;
  while((data[icnt] == ' ' || data[icnt] == '\n' || data[icnt] == '\t' || data[icnt] == '\v'
	 || data[icnt] == '\r' || data[icnt] == '\f') && icnt < MAXLINE)
    icnt++;
  for(i = 0; i < MAXLINE - icnt; i++)
    data[i] = data[i + icnt];

  /* checks that the line has not ended */
  if(data[0] == '\0' || icnt == MAXLINE)
    printf("Bad read error in init_read_ddatum.");

  /* puts the datum into the datum aray */
  i = 0;
  while((data[i] == '0' || data[i] == '1' || data[i] == '2' || data[i] == '3' ||
	 data[i] == '4' || data[i] == '5' || data[i] == '6' || data[i] == '7' ||
	 data[i] == '8' || data[i] == '9' || data[i] == '.' || data[i] == '-' ||
	 data[i] == '+' || data[i] == 'e' || data[i] == 'E' || data[i] == 'g' ||
	 data[i] == 'G') && i < MAXLINE)
    {
      datum[i] = data[i];
      i++;
    }

  /* checks that the next data is white space */
  if(data[i] != ' ' && data[i] != '\n' && data[i] != '\t' && data[i] != '\v' &&
     data[i] != '\r' && data[i] != '\f')
    printf("Bad double in init_read_ddatum.");

  /* peels off the data */
  for(j = 0; j < MAXLINE - i; j++)
    data[j] = data[i + j];

  /* reads the value */
  sscanf(datum, "%lf", &value);

  /* returns the value */
  return (value);
}


int init_read_idatum(
  char *line,			/* the original line */
  char *data			/* the data */
)
{
  int i, j;			/* loop counter */
  int icnt;			/* white space */
  char datum[MAXLINE];		/* holds the datum */
  int ivalue;			/* holds the value */

  /* initializes the datum */
  for(i = 0; i < MAXLINE; i++)
    datum[i] = '\0';

  /* strips the initial white space */
  icnt = 0;
  while((data[icnt] == ' ' || data[icnt] == '\n' || data[icnt] == '\t' || data[icnt] == '\v'
	 || data[icnt] == '\r' || data[icnt] == '\f') && icnt < MAXLINE)
    icnt++;
  for(i = 0; i < MAXLINE - icnt; i++)
    data[i] = data[i + icnt];

  /* checks that the line has not ended */
  if(data[0] == '\0' || icnt == MAXLINE)
    printf("Bad read error in init_read_idatum.\n");

  /* puts the datum into the datum aray */
  i = 0;
  while((data[i] == '0' || data[i] == '1' || data[i] == '2' || data[i] == '3' ||
	 data[i] == '4' || data[i] == '5' || data[i] == '6' || data[i] == '7' ||
	 data[i] == '8' || data[i] == '9' || data[i] == '-') && i < MAXLINE)
    {
      datum[i] = data[i];
      i++;
    }

  /* checks that the next data is white space */
  if(data[i] != ' ' && data[i] != '\n' && data[i] != '\t' && data[i] != '\v' &&
     data[i] != '\r' && data[i] != '\f')
    printf("Bad data in init_read_idatum.");

  /* peels off the data */
  for(j = 0; j < MAXLINE - i; j++)
    data[j] = data[i + j];

  /* reads the value */
  sscanf(datum, "%d", &ivalue);

  /* returns the value */
  return (ivalue);
}
