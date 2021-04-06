#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include<unistd.h>
#include<sys/stat.h>
#include<sys/types.h>

#include"rdtsc.h"

#define ALIGN 64

//Load all the FASTA datas
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

//Split datas in FASTA format
//1 file / sequence and namefile=IdSequence
//+ 1 file with all description lines
void split(const char *t){

  unsigned long long pos=0;
  unsigned long long t_len=strlen(t);
  
  unsigned long long seq_name_len=0;   //IdSequence length
  char *seq_name=NULL;                 //IdSequence
  unsigned long long seq_start=0;      //Position where a sequence start
  unsigned long long seq_len=0;        //Sequence length

  //File descriptor to open the file for sequence
  FILE *fd_s=NULL;
  
  unsigned long long des_line=0;       //Position where a description line start
  //Open the description file
  FILE *fd_d=fopen("description.txt","w");

  //Traversal all the datas, char by char
  for(pos=0;pos<t_len;pos++){
    //Check if we have a new sequence
    if(t[pos]=='>'){
      //Check if it's the first sequence
      if(seq_len>0){
	//Open sequence file
	fd_s=fopen(seq_name,"w");      
	//Store 1 sequence in her file
	fwrite(t+seq_start,sizeof(char),seq_len-2,fd_s);
	//Close the sequence file
	fclose(fd_s);
	//Free memory allocate for the IdSeq
	free(seq_name);
	fd_s=NULL;
      }
      //Reset sequence length
      seq_len=0;
      //Stock the description line start
      des_line=pos;
      //Go after '>'
      pos++;
      //Stock the sequence name start
      seq_name_len=pos;

      //Go to the first space
      while(t[pos]!=' ')
	pos++;

      //Allocate memory and stock IdSequence
      seq_name=strndup(t+seq_name_len,pos-seq_name_len);

      //Go to the next line
      while(t[pos]!='\n')
	pos++;

      //Store the description line
      fwrite(t+des_line,sizeof(char),pos-des_line+1,fd_d);

      //Stock the sequence start
      seq_start=pos+1;
    }

    //Count sequence length 
    seq_len++;
  }

  
  //Open latest sequence file
  fd_s=fopen(seq_name,"w");      
  //Store latest sequence in her file
  fwrite(t+seq_start,sizeof(char),seq_len,fd_s);
  //Close latest sequence file
  fclose(fd_s);
  //Free memory allocate for the IdSeq
  free(seq_name);
  
  //Close description file
  fclose(fd_d);
  
}


int main(int argc, char **argv){
  //Check arg
  if(argc<2)
    return printf("Usage: %s [file]\n",argv[0]);

  char *data=load_data(argv[1]);

  //Timestamp counter before the work
  double before=rdtsc();

  split(data);

  //Timestamp counter after the work
  double after=rdtsc();
  
  //Print the number of cycles to split
  printf("%lf cycles to split datas\n",after-before);

  release_data(data);

  return 0;
}
