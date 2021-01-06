#include <stdlib.h>
#include <stdio.h> 
#include <math.h> 
#include <time.h> 
#include <string.h>
#include <unistd.h>

#define NROWS 600
#define NCOLS 680
#define NSITES (NROWS*NCOLS)
#define TIMESTEPS (500)
#define PREYDIFF 3
#define LEFTBIAS (0.05)
#define TSKIP 10
#define XSKIP 1
#define esp 2

#define alpha 0.5
#define beta 1.0
#define epsilon 0.5
#define delta 0.1

//char* colors={"(0 \"red\", 1 \"green\")"}; 
char* colors={"(0 \"black\", 1 \"red\", 2 \"green\", 3 \"blue\", 4 \"yellow\")"};
FILE *gnuplot_pipe;
int preydiff;
double leftbias;
int **Mcell ;
int preycnt = 1400;
int predcnt = 800;
int itype = 0;
float NZERO = 0.3;
float AZERO = 0.2;
float wavelength = 10 ;


typedef struct siteinfo /* SITE INFO STRUCTURE DEFINITION */ {
  int prey, vbabe;
  int pred, pbabe;
  int indirect;
  int flag;
  struct siteinfo *nabes[4];
  } Site; // alias of struct siteinfo = Site 
Site habitat[NROWS][NCOLS];

void InitSeed(void){
	int seed;
	time_t nowtime; /* The time_t datatype is a data type in the ISO C library */
	struct tm *preztime; /* struct tm. Time structure */

	time(&nowtime);
	preztime = localtime(&nowtime);
	seed = (int)((preztime->tm_sec+1)*(preztime->tm_min+1)*(preztime->tm_hour+1)*(preztime->tm_year)*(preztime->tm_year));
	if(seed%2==0) seed++;
	srand48(seed);
}

void Initialize_Arena(Site *arena) {
 int i, j; 
 Site *origin;
 origin = arena; 
 for(i=0;i<NROWS;i++){
    for(j=0; j<NCOLS; j++){
        if(j>0) arena->nabes[0] = arena - 1; // arena->nabes[] acceder al vector nabes de la estructura tipo site
           else arena->nabes[0] = arena + NCOLS - 1;
        if(i>0) arena->nabes[1] = arena - NCOLS; 
           else arena->nabes[1] = arena + NSITES-NCOLS;
        if(i<(NROWS-1)) arena->nabes[2] = arena + NCOLS; 
           else arena->nabes[2] = origin + j;
        if(j<(NCOLS-1)) arena->nabes[3] = arena + 1; 
           else arena->nabes[3] = arena - (NCOLS-1);
        arena++;
    }
 }
}

void Initial_Conditions(Site arena[NROWS][NCOLS], int *nprey, int *npred){
	int i,j, total;
  Site *pnter;
  double cosval, hprob, pprob;
  pnter = &arena[0][0]; 

  for(i=0; i<NSITES; i++, pnter++)
          pnter->prey = pnter->pred = pnter->flag = 0;

  if(itype==0){
//      total= NZERO*NSITES; //NZERO es la poblacion inicial (densidad)
//      total = 4;
//      (*nprey) = *npred = total;
      pnter = &arena[0][0];
      for(i=0; i<preycnt; i++){
           j = drand48()*NSITES;
           while(pnter[j].prey==1) j = drand48()*NSITES;
      }     pnter[j].prey = 1;
      for(i=0; i<predcnt; i++){
           j = drand48()*NSITES;
           while(pnter[j].pred==1) j = drand48()*NSITES;
           pnter[j].pred = 1;
      }   
  }

   if(itype ==1){
       *nprey = *npred = 0;
       for(j=0; j<NCOLS; j++){
           cosval = cos(2.0*M_PI*j/wavelength);
           hprob = NZERO*(1.0+AZERO*cosval);
           pprob = NZERO*(1.0-AZERO*cosval);
           for(i=0;i<NROWS; i++){
                  if(drand48()<hprob){
                        arena[i][j].prey = 1;
                        (*nprey)++;
                  }
                  if(drand48()<pprob){
                        arena[i][j].pred = 1;
                        (*npred)++;
                  }
           }
       }
   }      
}

void DiffusePrey(Site *arena, double leftbias) {
	int i, dir;
	Site *origin, *nnabe;
	double xxx, bias[4];

	bias[0] = 0.25 + leftbias; bias[1] = bias [2] = 0.25;
	bias[3] = 0.25 - leftbias;

	for(i=0,origin = arena;i<NSITES;i++,origin++) {
        if(origin->prey==1){ /* consider move if present */
           xxx = drand48()-bias[0]; /* choose step direction */ 
           dir = 0; while(xxx>0.0) xxx -= bias[++dir];
           nnabe = origin->nabes[dir];

              if(nnabe->flag==0 && nnabe->prey==0) {/* empty and not tried before */
                  nnabe->indirect = 3 -dir;
                  nnabe->flag =1;  //flags step attempt 
                  origin->flag =2; /*individual moved */
              }    
              else if(nnabe->flag==1) {/* already stepped on */
                   nnabe ->flag = 3; /* flag mult steps */
                   nnabe ->nabes[nnabe->indirect]->flag = 0;
                   /* reset prev stepper */
              }     
        }
    }    
    origin = arena;
    for (i=0; i<NSITES; i++){
    	    if(origin->flag>0){
    	          if(origin->flag==1) origin->prey=1;
    	          else if (origin->flag==2)  origin->prey = 0;
    	        origin->flag = 0;
    	    }  
    	    origin++;
    }	         
}


void DiffusePred(Site *arena, double leftbias) {
  int i, dir;
  Site *origin, *nnabe;
  double xxx, bias[4];

  bias[0] = 0.25 + leftbias; bias[1] = bias [2] = 0.25;
  bias[3] = 0.25 - leftbias;

  for(i=0,origin = arena;i<NSITES;i++,origin++) {
        if(origin->pred==1){ /* consider move if present */
           xxx = drand48()-bias[0]; /* choose step direction */ 
           dir = 0; while(xxx>0.0) xxx -= bias[++dir];
           nnabe = origin->nabes[dir];

              if(nnabe->flag==0 && nnabe->pred==0) {/* empty and not tried before */
                  nnabe->indirect = 3 -dir;
                  nnabe->flag =1;  //flags step attempt 
                  origin->flag =2; /*individual moved */
              }    
              else if(nnabe->flag==1) {/* already stepped on */
                   nnabe ->flag = 3; /* flag mult steps */
                   nnabe ->nabes[nnabe-> indirect]->flag = 0;
                   /* reset prev stepper */
              }     
        }
    }    
    origin = arena;
    for (i=0; i<NSITES; i++){
          if(origin->flag>0){
                if(origin->flag==1) origin->pred=1;
                else if (origin->flag==2)  origin->pred = 0;
              origin->flag = 0;
          }  
          origin++;
    }          
}


void Interact(Site *arena, int *nprey, int *npred){
	int i, dir, flag, bcnt;
	Site *origin, *nnabe;
	origin = arena; /* prey production */
	for(i=0; i<NSITES; i++, origin++){
		if(origin->prey==1){
			  if(drand48()<alpha){
			  	    nnabe = origin;
			  	    flag = bcnt = 1;
			  	    while(flag==1 && bcnt<=10){
			  	    	dir = 4*drand48();
			  	    	nnabe = nnabe->nabes[dir];
			  	    	if (nnabe->prey==1 || nnabe->vbabe ==1){
			  	    		bcnt++;
			  	    	}
                        else {
                        	nnabe->vbabe = 1;
                        	flag = 0;
                        }
			  	    }
			  }
		}
	}
	/*predation and predator reproduction */
	origin = arena; 
	for(i=0; i<NSITES; i++, origin++){
		      if(origin->prey * origin->pred==1){
		      	          origin->prey =0;
		      	          if (drand48()<epsilon){
		      	          	         nnabe = origin;
		      	          	         flag = bcnt = 1;
		      	          	         while(flag ==1 && bcnt <=10){
		      	          	         	       dir = 4*drand48();
		      	          	         	       nnabe = nnabe->nabes[dir];
		      	          	         	       if ( nnabe->pred==1 || nnabe->pbabe==1){
		      	          	         	       	            bcnt++;
		      	          	         	       }
		      	          	         	       else {
		      	          	         	       	nnabe->pbabe =1 ;
		      	          	         	       	flag = 0;

		      	          	         	       }

		      	          	         }
		      	          }
		      }
	}

	origin = arena;  /* predator death */
	for(i=0; i>NSITES; i++, origin ++){
		     if(origin->pred==1){
		     	      if (drand48()<delta){
		     	      	origin->pred = 0;
		     	      }

		     }
	}

	*nprey = *npred = 0;
	origin = arena;
	for(i=0; i<NSITES; i++, origin++){
		     if(origin->vbabe==1){
		     	       origin->vbabe = 0;
		     	       origin ->prey = 1;	       
		     } 
		     if(origin->pbabe==1){
		     	       origin->pbabe = 0;
		     	       origin->pred = 1;
		     }
		     if(origin->prey) (*nprey)++;
		     if(origin->pred) (*npred)++;

	}

}


void fill_matrix_pointer(Site *arena){
  int i, j =0;
  int cero =0;
  int uno = 1;
  Site *origin;
  origin = arena ;
   
     for (i=0; i<NROWS; i++){
        for(j=0; j<NCOLS; j++){
           if(arena->prey==0 && arena->pred==0)
                Mcell[i][j] = 0;
            else if(arena->prey==1)
                  Mcell[i][j] = 1;
            else if(arena->pred==1) 
                   Mcell[i][j] = 2;
            arena++;     
        }
//        printf("\n"); 
     }
}

int main(void) {
    int i, j, step;
    Site *plataform;
    FILE *imfileid;
    preydiff = PREYDIFF;
    leftbias = LEFTBIAS;

    Mcell=(int**) calloc(NROWS+2,sizeof(int*));
    for(i=0;i<NROWS+2;i++)
    Mcell[i]=(int*) calloc(NCOLS+2,sizeof(int));
  
/// pipe para configurar gnuplot
    gnuplot_pipe = popen("gnuplot", "w");
    fprintf(gnuplot_pipe, "set xrange [0:%d]\n",NROWS);
    fprintf(gnuplot_pipe,"set yrange [0:%d] reverse\n",NCOLS);
    fprintf(gnuplot_pipe,"set size ratio 1.0\n");
    fprintf(gnuplot_pipe,"set cbrange [0:%d]\n",esp);
    fprintf(gnuplot_pipe,"set palette defined %s\n",colors);
    fprintf(gnuplot_pipe,"unset colorbox\n");
  //fprintf(gnuplot_pipe,"set pm3d map\n");
  //fprintf(gnuplot_pipe,"set pm3d interpolate 2,2\n"); */

    InitSeed();
    Initialize_Arena(&habitat[0][0]);
    Initial_Conditions(habitat,&preycnt, &predcnt);
 //   my_initial_conditions(habitat, preycnt, preydiff);
 //   other_initial_conditions(habitat, preycnt, preydiff);


       for(step=0;step<TIMESTEPS;step++) {

         	for(i=0;i<preydiff;i++)
             DiffusePrey(&habitat[0][0], leftbias);
             DiffusePred(&habitat[0][0], leftbias); 
             Interact(&habitat[0][0], &preycnt, &predcnt);
 

         //   fprintf(gnuplot_pipe,"set title \"k = %d MCS\"\n",step);
         //   fprintf(gnuplot_pipe,"plot '-' matrix with image t ''\n");
  

     //            online_print_matrix(habitat);

      //      print_matrix(habitat, gnuplot_pipe);

       //        online_transform_matrix(habitat);

      //       transform_matrix(habitat, gnuplot_pipe);
               
   //            fprintf(gnuplot_pipe,"\n"); 

   //         fill_matrix(habitat);

   //         fill_two_matrix(habitat);

  //           print_sites(habitat);

            fill_matrix_pointer(&habitat[0][0]);

            fprintf(gnuplot_pipe,"set title \"k = %d Monte Carlo Steps\"\n",step);
            fprintf(gnuplot_pipe,"plot '-' matrix with image t ''\n");  
            for (i=0;i<=NROWS+1;i++) {
                for (j=0;j<=NCOLS+1;j++){
                    fprintf(gnuplot_pipe,"%d ", Mcell[i][j]);
                    //printf("%d ", Mcell[i][j]);
                }
        //        printf("\n");
                fprintf(gnuplot_pipe,"\n");
            }
            fprintf(gnuplot_pipe,"\n e");
            fprintf(gnuplot_pipe,"pause 0.5\n");  
  
     

       //         fprintf(gnuplot_pipe,"\n e");
       //         fprintf(gnuplot_pipe,"pause 0.5\n");    
             }  
           fclose(gnuplot_pipe);  
          // fclose(imfileid);
 
   return 0;
}


