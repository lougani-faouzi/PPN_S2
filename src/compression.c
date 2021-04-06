#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define CODE_A 0
#define CODE_T 3
#define CODE_C 1
#define CODE_G 2

//
int main(int argc, char **argv){
  //Check arg
  if(argc<2)
    return printf("Usage: %s file_sequence\n",argv[0]), 1;

  //Open files, compressed and not
  FILE *fd=fopen(argv[1], "r");

  //Create the filename of compressed file
  char *filename=strcat(argv[1],".compressed");
  
  FILE *fd_c=fopen(filename, "w");

  //Check pointers
  if(!fd)
    return printf("Error: can't open file sequence\n"), 2;
  if(!fd_c)
    return printf("Error: can't open file compressed\n"), 2;

  unsigned char compress=0;
  int i=0;
  //Traversal file sequence char by char
  char c;
  while((c=fgetc(fd))!=EOF){
    //Compress 4 octets in 1 octet(8bits in 2bits)
    switch(c){
    case 'A':
      compress=compress<<2;
      compress=compress | CODE_A;
      break;
    case 'T':
      compress=compress<<2;
      compress=compress | CODE_T;
      break;
    case 'C':
      compress=compress<<2;
      compress=compress | CODE_C;
      break;
    case 'G':
      compress=compress<<2;
      compress=compress | CODE_G;
      break;
    default:
      continue;
    }
    
    i++;
    if(i>3){
      fputc(compress,fd_c);
      compress=0;
      i=0;
    }
  }

  while(i++<4)
    compress=compress<<2;
  fputc(compress,fd_c);
  
  //Close all files
  fclose(fd);
  fclose(fd_c);
  
  return 0;
}
