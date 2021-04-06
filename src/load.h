#include<sys/stat.h>
#include<sys/types.h>
#include<stdio.h>
#include<stdlib.h>

#define ALIGN 64

//Load all the data
char *load_data(const char *fname){
  //check arg
  if(!fname)
    return printf("Error: null pointer!\n"),NULL;

  //Structure to store meta-info
  struct stat sb;

  //Load meta-info
  if(stat(fname,&sb)<0)
    return printf("Error: can't stat %s !\n",fname),NULL;

  //Open file data
  FILE *fd=fopen(fname,"r");

  //check pointer
  if(!fd)
    return printf("Error: can't open %s !\n",fname),NULL;

  //Allocate memory for store datas
  char *t= aligned_alloc(ALIGN,sb.st_size+1);

  //Check pointer
  if(!t)
    return printf("Error: can't allocate memory !\n"),NULL;

  //Load all datas
  unsigned long long n= fread(t,sizeof(char),sb.st_size,fd);

  //Close file data
  fclose(fd);

  //Check no missing data
  if(n != sb.st_size)
    return printf("Error: read size doesn't match !"),NULL;

  //Make sure data is zero terminated
  t[sb.st_size]=0;

  return t;
}

//Free datas memory
static inline void release_data(char *t){
  //Check arg
  if(t)
    free(t);
}


