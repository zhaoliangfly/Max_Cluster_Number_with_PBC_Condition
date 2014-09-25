#include <stdio.h>
#include <stdlib.h>
#include "xdrfile.h"
#include "math.h"
#include "time.h"
#include "Cluster_Parameters_SU5_OH.h"
#include "Cluster_Functions_SU5_OH.h"

#define Specify_Cluster_Max_Number_Search_CM   0
#define Max_Cluster_Number      1
#define All_Cluster_Number      0



#define TIME 0
#define TIME_FI  200000
#define Deltar   0.001
#define Dotsnumber  2200 
#define Pdf_Cutoff   2.2 
#define Cluster_Atoms_Distance   0.60
#define Cluster_CM_SU5_OH_Distance  0.9
#define Boundary_Margin 0.7

void SEAL(int i);

int natoms,step,natom,read_return;
float time_1,prec;
//float prec;
matrix box;
rvec *x;
XDRFILE *xtc;
 
int cluster;
int i,j;
          struct data
          {
                 int label;
                 int location_x;
                 int location_y;
                 int location_z;
          };
         struct data array[Molnumber_SU5_OH];   /** Notice in this array, the array[0] is also used to store the data! **/


int main ()
{
   long t1;
   t1=clock();
   
   char path_xtc[]="F:\\GMXSimulationData\\SuperComputation6202WatersRealSU5\\42\\traj.xtc";

    xtc=xdrfile_open (path_xtc,"r");
    read_xtc_natoms (path_xtc,&natoms);
    x = calloc(natoms, sizeof (x[0]));


  int time_initial=0;  /* The starting point of our investigation depends on the real time of the .xtc file! **/
  int time_statistic=TIME;
  double pdistribution[2201];
    for(i=0;i<=Dotsnumber;i++)
          pdistribution[i]=0.00;
  
  double All_Cluster[Molnumber_SU5_OH+1];
    for(i=0;i<=Molnumber_SU5_OH;i++)
         All_Cluster[i]=0.00;


while (1)
{
  read_return=read_xtc (xtc,natoms,&step,&time_1,box,x,&prec);
            if (read_return!=0) break;
   

  
#if 0
    for (natom=1;natom<=natoms;natom++)
    {
    printf ("%d %f %d %f %f %f\n",step,time,natom,x[natom-1][0],x[natom-1][1],x[natom-1][2]);
    }
#endif  
  
   

   
   if(time_statistic==time_initial)  
   {
         j=-1;
         for(i=0;i<=Molnumber_SU5_OH*Atomnumber_SU5_OH-1;i+=1)
         {
              if(i%Atomnumber_SU5_OH==0)
                   j++;
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].x=x[i][0];
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].y=x[i][1];
              MOL_SU5_OH_48[j].MOL_SU5_OH_1[i%Atomnumber_SU5_OH].z=x[i][2];
         }

         j=-1;
         for(i=Molnumber_SU5_OH*Atomnumber_SU5_OH;i<=Molnumber_SU5_OH*Atomnumber_SU5_OH+3*Molnumber_SOL-1;i+=1)
         {
              if((i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3==0)
                   j++;
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].x=x[i][0];
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].y=x[i][1];
              MOL_SOL_3559[j].MOL_SOL_1[(i-Molnumber_SU5_OH*Atomnumber_SU5_OH)%3].z=x[i][2];
          }



         for(i=0;i<=Molnumber_SU5_OH-1;i++)
         {
               array[i].label=0;
               array[i].location_x=0;
               array[i].location_y=0;
               array[i].location_z=0;
         }
     /** For the convienience of Computation ***
      ** we define the CM matrix explicitly  ***
      **/
         Position CM_Of_SU5_OH[Molnumber_SU5_OH];
           for(i=0;i<=Molnumber_SU5_OH-1;i++)
           {
                 CM_Of_SU5_OH[i].x=CM_SU5_OH_X(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
                 CM_Of_SU5_OH[i].y=CM_SU5_OH_Y(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
                 CM_Of_SU5_OH[i].z=CM_SU5_OH_Z(&(MOL_SU5_OH_48[i].MOL_SU5_OH_1[0]));
           }

           
         
         Position CM_Of_SOL[Molnumber_SOL];
           for(i=0;i<=Molnumber_SOL-1;i++)
           {
               CM_Of_SOL[i].x=CM_SOL_X(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
               CM_Of_SOL[i].y=CM_SOL_Y(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
               CM_Of_SOL[i].z=CM_SOL_Z(&(MOL_SOL_3559[i].MOL_SOL_1[0]));
            }
         
         


             cluster=0;
              for(j=0;j<=Molnumber_SU5_OH-1;j++)
              {
                   if(array[j].label==0)   /***To prevent the repeated labeling******/
                        {
                           cluster++;
                           SEAL(j); 
                         }  
              }

    /** End of the Cluster Labeling **/
      
    #if 0
  
      for(i=0;i<=Molnumber_SU5_OH-1;i++)
               printf("    1SU5  A   %5d%8.3f%8.3f%8.3f\n",i+1,CM_Of_SU5_OH[i].x,CM_Of_SU5_OH[i].y,CM_Of_SU5_OH[i].z);

      printf("   %8.3f   %8.3f   %8.3f\n", box[0][0], box[1][1], box[2][2]);

    #endif

#if  Max_Cluster_Number

    int Recording_Cluster[Molnumber_SU5_OH+1];
    int Recording_Label;    /** This number can tell us which label corresponds to the max Cluster Number! **/
    for(i=0;i<=Molnumber_SU5_OH;i++)
        Recording_Cluster[i]=0;

    for(i=0;i<=Molnumber_SU5_OH-1;i++)
        Recording_Cluster[array[i].label]++;
    int max=1;
    for(i=0;i<=Molnumber_SU5_OH;i++)
    {
          if(Recording_Cluster[i]>=max)
          {
                max=Recording_Cluster[i];
                Recording_Label=i;
           }
     }

     printf("%d %d\n",time_statistic, max);
#endif

#if  All_Cluster_Number

    int Recording_Cluster[Molnumber_SU5_OH+1]; /** The label begins from 1 ! **/
    int Recording_Label;    /** This number can tell us which label corresponds to the max Cluster Number! **/
    for(i=0;i<=Molnumber_SU5_OH;i++)
        Recording_Cluster[i]=0;

    for(i=0;i<=Molnumber_SU5_OH-1;i++)
        Recording_Cluster[array[i].label]++;

    for(i=0;i<=Molnumber_SU5_OH;i++)
    { 
    	if(Recording_Cluster[i]!=0)	
        All_Cluster[Recording_Cluster[i]]++;
    }    
      

        
#endif






   
      

           /** Here is the end of the labeling! **/

           time_statistic++;
          if(time_statistic==TIME_FI+1)
               break;

    }  /** End of the time_statistic ***/

    time_initial++;

}   /** End of the loops ***/

#if All_Cluster_Number
      for(i=0;i<=Molnumber_SU5_OH;i++)
      {
      	 if(All_Cluster[i]!=0)
           printf("%d %f\n",i,All_Cluster[i]/((double)(TIME_FI-TIME)));
      }     
           
#endif           
  
    xdrfile_close (xtc);
    free(x);

    //getchar();

    long t2;
    t2=clock();
    printf("%d\n",t2-t1);

        return 0;


}   /** End of the main() ***/

     /****The sub-code below is the coral codes for labeling********/    
         void SEAL( int s)    
         {
             array[s].label=cluster;
             int neighbour[Molnumber_SU5_OH+1];
             for(i=0;i<=Molnumber_SU5_OH;i++)
                neighbour[i]=-1;
             int u=0;

             for(i=0;i<=Molnumber_SU5_OH-1;i++)    /****We try to judge whether  the molecules can be relative to objetive atom*****/ 
             {
                   if(i==s)
                       continue;
                   else
                   {
                       int px;
                       int py;
                       int pz;
                       int PX;
                       int PY;
                       int PZ;
                       double Dis=20.0;

                       int v1,v2;
                   for(v1=0;v1<=3;v1++)
                   {
                     for(v2=0;v2<=3;v2++)
                     {

                       for(px=-1;px<=1;px++)    /****We must take the peoredical boxes into considerations******/
                       {
                          for(py=-1;py<=1;py++)
                          {
                             for(pz=-1;pz<=1;pz++)
                             {

                                   if(distance(MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].x,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].y,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].z,MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].x+px*box[0][0],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].y+py*box[1][1],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].z+pz*box[2][2])<=Dis)
                                   {
                                           Dis=distance(MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].x,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].y,MOL_SU5_OH_48[s].MOL_SU5_OH_1[v1].z,MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].x+px*box[0][0],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].y+py*box[1][1],MOL_SU5_OH_48[i].MOL_SU5_OH_1[v2].z+pz*box[2][2]);
                                   }

                               }
                            }
                         }
                        }
                       }
       
                       
                                         
                                   if((Dis<=Cluster_Atoms_Distance)&&(array[i].label==0))
                   /****The condition at the last can ensure the finity of loops. And it is very important!**************/                  
                                   {
                                          neighbour[u]=i;
                                          array[i].label=cluster;
                                          u++;
                                          
                                   }
           
                        }
                     } 

         
               int k;
               
               for(k=0;k<=Molnumber_SU5_OH;k++)
               {      
                      if(array[neighbour[k]].label==-1)
                               break;
                      if(array[neighbour[k]].label==cluster)
                               SEAL(neighbour[k]);
               }
                
          }
