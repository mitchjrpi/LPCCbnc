#include <ilcplex/cplex.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "LPCC.h"
#include "LPCCsolver.h"
#include <time.h>
#include "rngs.h"
#include "rvgs.h"
#include "utilities.h"
//#include <vld.h>
int nodecnt = 0;
int leftnodes = 0;
int main (int argc, 
		  char **argv)
{
   int status = 0;
   PARAM inputParam;
   parse_command_line(argc,argv,&inputParam);
   status=solveLPCC(inputParam);
   if(status) goto TERMINATE;
TERMINATE:
   return (status);

}  

