#define MAX_GENES 1000

typedef struct gene_map_s {

   //
   unsigned long long gene_counter;

   //Gene start position (ATG)
   unsigned long long gene_start[MAX_GENES];

   //Gene stop position (TAA, TAG, TGA)
   unsigned long long gene_end[MAX_GENES];
}* gene_map;


int min(int a, int b){
  if(a<=b)
    return a;
  else
    return b;
}


int find(char* motif,char* seq){
  int boolean=0;
  int pos=0;
  while(!boolean && seq[pos+2]!='\0'){
    if(motif[0]==seq[pos]){
      if(motif[1]==seq[pos+1]){
        if(motif[2]==seq[pos+2]){
          boolean=1;
        }
      }
    }
    pos+=3;
  }
  if(boolean)
    return pos;
  return -1;
}

int find_stop(char* seq){
  int stop1=0;
  int stop2=0;
  int stop3=0;

  stop1=find("TAA",seq);
  stop2=find("TAG",seq);
  stop3=find("TGA",seq);

  //Si je trouve pas TAA
  if(stop1<0){
    //Si je trouve pas TAA ni TAG
    if(stop2<0){
      //Je retourne la pos TGA (-1 si rien trouvé, sa pos sinon)
      return stop3;
    //Si je trouve pas TAA mais je trouve TAG
    }else{
      //Si je trouve pas TGA
      if(stop3<0){
        //Je n'ai trouvé que TAG donc je retourne sa pos
        return stop2;
      }else{
        //J'ai trouvé TAG et TGA donc je choisis le minimum des deux
        return min(stop2,stop3);
      }
    }
  //Si je trouve TAA
  }else{
    //Si je trouve pas TAG
    if(stop2<0){
      //Si je trouve pas TGA
      if(stop3<0){
        //Je n'ai trouvé que TAA je retourne sa pos
        return stop1;
      //Si je trouve TGA
      }else{
        //J'ai trouvé TAA et TGA donc je retourne la pos minimale
        return min(stop1,stop3);
      }
    //Si je trouve TAG
    }else{
      //Si je trouve pas TGA
      if(stop3<0){
        //J'ai trouvé TAA et TAG donc je retourne la pos minimale
        return min(stop1,stop2);
      //Si je trouve TGA
      }else{
	//Je trouve TAA, TAG et TGA, je retourne le min entre les 3
        int mini=min(stop1,stop2);
        return min(mini,stop3);
      }
    }
  }
}


