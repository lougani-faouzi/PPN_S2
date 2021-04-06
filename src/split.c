#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include"rdtsc.h"

//Number of character max in 1 line
//Can be found with "cat sequences.fasta | grep "^>" | wc -L
#define BUFFER_SIZE 600

#define KEY_LENGTH 10

int main(int argc, char **argv){
  //Check arg
  if(argc<2)
    return printf("Usage: %s file_data\n",argv[0]);

  //Open file data
  FILE *fd=fopen(argv[1],"r");
  //Create file to store sample description
  FILE *fd_d=fopen("description.txt","w");
  //Check pointers
  if(!fd)
    return printf("Error: can't open file data\n"), 1;
  if(!fd_d)
    return printf("Error: can't open file description\n"), 1;


  //Take timestamps counter before work
  double before=rdtsc();
  
  //Allocate memory for buffer
  char line[BUFFER_SIZE];
  
  //File descriptor to store sequences
  FILE *fd_s=NULL;
  
  int pos;
  char *file_name=NULL;
  
  //Traversal file data line by line
  while(fgets(line,BUFFER_SIZE,fd)!=NULL){
    if(line[0]=='>'){

      fputs(line,fd_d);
      fputc('\n',fd_d);
      
      //Check if a file sequence is open
      if(fd_s){
	fclose(fd_s);
	free(file_name);
      }

      pos=0;
      while(line[pos]!=' ')
	pos++;
      
      file_name=strndup(line+1,pos-1);

      //Open a file sequence
      fd_s=fopen(file_name,"w");
      if(!fd_s)
	return printf("Error: can't open files sequence"),1;

    }else{
      fputs(line,fd_s);
    }
  }

  //Take timestamps counter after work
  double after=rdtsc();

  //Calcul the time to split
  double time=after-before;

  //Print the time to split the datas
  printf("%lf cycles to split datas\n",time);


  //Close all files
  fclose(fd);
  fclose(fd_s);
  fclose(fd_d);
  return 0;
}
