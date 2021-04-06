#include<unistd.h>
#include<string.h>

#include"load.h"
#include"detection.h"
#include"popcount.h"
#include"rdtsc.h"


#define SEUIL 10


//Genere la mapping des codons start et stop d'une sequence ADN
gene_map mapping(char *seq){

  //Alloue la memoire pour la map des genes
  gene_map gm=malloc(sizeof(struct gene_map_s));
  //Initialise le compteur de gene a 0
  gm->gene_counter = 0;

  //Variable pour enregistrer la pos d'un codon start
  //Seulement si y en a un de trouver
  int start_pos = 0;
  //Variable pour enregistrer la pos d'un codon stop
  //Seulement si y en a un de trouver
  int stop_pos = 0;

  //Variable pour check si on trouve notre motif
  int check=-1;

  //Chercher les zones codantes / mapper la sequence
  do{
    //Cherche le motif "ATG" dans la sequence à partir de la position du dernier codon stop trouvé
    //stop_pos: 0 à l'init
    check=find("ATG",seq+stop_pos);
    //Si on a trouvé le motif
    if(check!=-1){
      //On enregistre sa position dans la sequence
      start_pos=stop_pos+check-3;
      /*printf("start_pos : %d "
               "-> codon : %c%c%c\n",start_pos,seq[start_pos],seq[start_pos+1],seq[start_pos+2]);
      */

      //Cherche le premier codon stop à partir de la position du dernier codon start trouvé
      check=find_stop(seq+start_pos);
      //Si on a trouvé un stop
      if(check!=-1){
        //On enregistre sa position dans la sequence
        stop_pos=start_pos+check-3;
        /*printf("stop_pos : %d "
                 "-> codon : %c%c%c\n",stop_pos,seq[stop_pos],seq[stop_pos+1],seq[stop_pos+2]);
        */

        //On enregistre les positions trouvé dans notre structure, map
        gm->gene_start[gm->gene_counter]=start_pos;
        gm->gene_end[gm->gene_counter]=stop_pos+2;
        gm->gene_counter++;
        //printf("ADD !\n");
      }
    }
  //Tant que l'on trouve un codon start et stop
  }while(check!=-1);

  return gm;
}


//Transnforme l'ADN en ARN messager, seulement pour les genes de la séquence
char **generate_ARN(gene_map gm,char* seq){

  //Alloue de la memoire pour stocker les ARN dans un tableau
  char** ARN_m = (char **)malloc(gm->gene_counter*sizeof(char*));
  //Pour chaque ligne on alloue l'espace nécessaire pour stocker le gene en ARN
  for(int i=0; i < gm->gene_counter; i++)
    ARN_m[i] = (char *)malloc(gm->gene_end[i] - gm->gene_start[i] + 1 *sizeof(char));


  int pos=0;
  int pos_arn=0;
  //Pour chaque gene trouvé
  for(int i=0; i < gm->gene_counter; i++){
    //Init la pos au debut du gene
    pos=gm->gene_start[i];

    pos_arn=0;

    //On parcours le gene caractére par caractére pour l'afficher sous forme ADN
    while(pos<=gm->gene_end[i]){
      //Si on croise le caractère T on le remplace par U
      if(seq[pos]=='T')
        ARN_m[i][pos_arn]='U';
      else
        ARN_m[i][pos_arn]=(char)seq[pos];

      pos++;
      pos_arn++;
    }
  }

  return ARN_m;
}

//
void generate_prot(char *arn, char **codons){
  int i=0;
  int pos=0;
  char prot[3];
  //Nombre proteines
  int nb_prot=strlen(arn)/3;
  int short_names[nb_prot];
  int symbols[nb_prot];

  printf("ARN_m : ");
  
  //Parcours chaque codons de la sequence ARN
  for(pos=0; pos<nb_prot; pos++){
  	//Recupere le codon dans prot
    prot[0]=arn[pos*3];
    prot[1]=arn[pos*3+1];
    prot[2]=arn[pos*3+2];

    //Compare a tout les codons possible
    //pour trouver a quelle proteine il correspond
    i=0;
    while(strcmp(codons[i],prot) && i<192){
      i+=3;
    }
    //Une fois trouver on enregistre la position du short name
    short_names[pos]=i+1;
    //Une fois trouver on enregistre la position du symbol
    symbols[pos]=i+2;
    printf("%s ",prot);
  }
  //Affichage en parcourant chaque proteine obtenu
  printf("\n(Short name) Protein : ");
  for(pos=0; pos<nb_prot; pos++)
    printf("%s ",codons[short_names[pos]]);
  printf("\n(Symbol) Protein : ");
  for(pos=0; pos<nb_prot; pos++)
    printf("%s ",codons[symbols[pos]]);
  printf("\n");

}

//
void detect_mut(char *gene){
  //Position dans le gene
  int pos=0;
  //Le nombre de G et de C a la suite
  int occurrence=0;
  //Le % de G et C à la suite par rapport aux nombres totals de nucléotide
  double risk_rate=0.00;
  //Nombre de zone a risque trouve
  int i=0;
  //Compteur de Zones (qui sont POTENTIELLEMENT a risque)
  int cpt_zar=0;
  //Tableau pour stocker lesquelles sont a risque parmis toute les zones POTENTIELLES
  int zar[1000]={0};

  //ON parcour le gene caractere par caractere
  while(gene[pos]!='\0'){

    //Si on croise un G ou un C, on incrémente l'occurrence
    if(gene[pos]=='G' || gene[pos]=='C'){
      if(occurrence==0)
        //Début d'une zone
        printf("[");

      occurrence++;

    //Sinon
    }else{
      //Si on a croisé une zone
      if(occurrence>0){
        //Fin d'une zone
        printf("]");
        //On incrémente le compteur de zone
        cpt_zar++;

        //On calcule le %
        risk_rate = occurrence * 100.00;
        risk_rate = risk_rate / strlen(gene);

        //Si on depasse le seuil pour etre a risque
        if(risk_rate>SEUIL){
          //On a trouvé une zone a risque qui commence en pos-occurence
          zar[i]=cpt_zar;
          //On incrémente le nombre de zone a risque
          i++;
        }
        
        //On réinitialise l'occurrence pour terminer la zone
        occurrence=0;
      }

    }
    //Affiche le caractere courant
    printf("%c",gene[pos]);
    //On passe au caractere suivant
    pos++;

  }
  printf("\n");

  //Si le 1er entier est a 0, on a trouve aucune zone a risque
  i=0;
  if(zar[0]==0)
    printf("\nAucune zone est à risque de mutation\n");
  //Sinon on les affiches
  while(zar[i]!=0)
    printf("\nz%d est à risque de mutation\n", zar[i++]);

  printf("\n");

}


//Calcul le pourcentage de correspondance entre les 2 sequences
//seq DOIT etre la plus grande
void matching_rate(char *seq, char *seq2){
  //On stocke le minimum entre les tailles des 2 sequences
  int mini=min(strlen(seq),strlen(seq2));

  //
  double d=0;
  int pos=0;
  //Parcours la plus grande sequence, de base en base
  while(seq[pos+mini-1]!='\0'){
    //Affiche le morceau de la grande sequence qui sera comparee 
    for(int i=0;i<mini;i++)
      printf("%c",seq[pos+i]);

    printf("\n");
    //Affiche la petite sequence qui sera comparee
    for(int i=0;i<mini;i++)
      printf("%c",seq2[i]);

    //XOR caractere par caractere les deux chaines
    //Puis compte le nombre de bits a 1 total 
    for(int i=0;i<mini;i++)
      d += popcount( seq[pos+i] ^ seq2[i] );

    //Nombre total de bits -> 100% de bits differents
    //Nombre de bit a 1    -> y% de bits differents
    d = d * 100;
    d = d / (8*mini);

    //Affiche le résultat trouve
    printf("\n\nIl y a %lf%% de bits diférrents entre ces 2 séquences.\n\n",d);

    d=0;
    pos++;
  }
}


//
int main(int argc, char **argv){
  //Check arg
  if(argc<3)
    return printf("Usage: %s [file seq1] [file seq2]\n",argv[0]);

//Charge les codons en mémoire
char* codons[192]={"AAA","Lys","K",
		    		 "AAC","Asn","N",
		    		 "AAG","Lys","K",
		     		 "AAU","Asn","N",
					 "ACA","Thr","T",
					 "ACC","Thr","T",
					 "ACG","Thr","T",
					 "ACU","Thr","T",
					 "AGA","Arg","R",
					 "AGC","Ser","S",
					 "AGG","Arg","R",
					 "AGU","Ser","S",
					 "AUA","Ile","I",
					 "AUC","Ile","I",
					 "AUG","Met","M",
					 "AUU","Ile","I",
					 "CAA","Gln","Q",
					 "CAC","His","H",
					 "CAG","Gln","Q",
					 "CAU","His","H",
					 "CCA","Pro","P",
					 "CCC","Pro","P",
					 "CCG","Pro","P",
					 "CCU","Pro","P",
					 "CGA","Arg","R",
					 "CGC","Arg","R",
					 "CGG","Arg","R",
					 "CGU","Arg","R",
					 "CUA","Leu","L",
					 "CUC","Leu","L",
					 "CUG","Leu","L",
					 "CUU","Leu","L",
					 "GAA","Glu","E",
					 "GAC","Asp","D",
					 "GAG","Glu","E",
					 "GAU","Asp","D",
					 "GCA","Ala","A",
					 "GCC","Ala","A",
					 "GCG","Ala","A",
					 "GCU","Ala","A",
					 "GGA","Gly","G",
					 "GGC","Gly","G",
					 "GGG","Gly","G",
					 "GGU","Gly","G",
					 "GUA","Val","V",
					 "GUC","Val","V",
					 "GUG","Val","V",
					 "GUU","Val","V",
					 "UAA","Stp","O",
					 "UAC","Tyr","Y",
					 "UAG","Stp","O",
					 "UAU","Tyr","Y",
					 "UCA","Ser","S",
					 "UCC","Ser","S",
					 "UCG","Ser","S",
					 "UCU","Ser","S",
					 "UGA","Stp","O",
					 "UGC","Cys","C",
					 "UGG","Trp","W",
					 "UGU","Cys","C",
					 "UUA","Leu","L",
					 "UUC","Phe","F",
					 "UUG","Leu","L",
					 "UUU","Phe","F"};


  printf("*** CHARGE LA PREMIERE SEQUENCE ***\n\n");
  //Charge la sequence
  char *seq=load_data(argv[1]);


/// MAPPING
  
  //
  double before=rdtsc();

  printf("*** MAPPING EN COURS ***\n\n");
  //Map la sequence avec la structure
  gene_map gm=mapping(seq);  
  
  //
  double after=rdtsc();
  printf("%lf cycles to find all genes\n\n",after-before);


/// GENERE ET AFFICHE TOUT LES GENES SOUS FORME D'ARN MESSAGER

  //
  before=rdtsc();

  printf("*** GENERATION DE L'ARN MESSAGER POUR CHAQUE GENE ***\n\n");
  //Genere l'ARN_m pour chaque gene de la sequence
  char** ARN_m=generate_ARN(gm,seq);

  //
  after=rdtsc();

  //
  int pos=0;
  //Pour chaque gene trouvé
  for(int i=0; i < gm->gene_counter; i++){
    printf("(ADN)Gene %d : ",i+1);

    //Initialise la pos au debut du gene
    pos=gm->gene_start[i];

    //On parcours le gene caractére par caractére pour les afficher
    while(pos<=gm->gene_end[i]){
      printf("%c", seq[pos]);
      pos++;
    }
    printf("\n");

    printf("(ARN)Gene %d : ",i+1); 
    printf("%s",ARN_m[i]);  
    printf("\n\n");
  }

  printf("%lf cycles to generate mARN for all genes\n\n",after-before);


/// GENERE ET AFFICHE LES PROTEINES SELON L'ARN
  //
  before=rdtsc();

  printf("*** GENERATION DES PROTEINES ***\n\n");
  //Genere l'ARN_m pour chaque gene de la sequence
  for(int i=0; i < gm->gene_counter; i++){
    generate_prot(ARN_m[i], codons);
    printf("\n");
  }

  //
  after=rdtsc();
  printf("%lf cycles to generate protein for all genes\n\n",after-before);


///DETECTION DES ZONES A RISQUES DE MUTATIONS

  printf("*** RECHERCHE DE ZONE A RISQUE DE MUTATION ***\n\n");
  //Lance la detection de zone a risque de mutation pour chaque gene de la sequence
  for(int i=0; i < gm->gene_counter; i++){
    detect_mut(ARN_m[i]);
    printf("\n");
  }
  printf("\n");
  //Libère la mémoire des genes en ARN
  free(ARN_m);
  //Libère la mémoire du mapping
  free(gm);


///MATCHING RATE SUR LES BITS (PAS SUR LES BASES)

  printf("*** CHARGE LA DEUXIEME SEQUENCE ***\n\n");
  //Charge la sequence
  char *seq2=load_data(argv[2]);

  matching_rate(seq,seq2);


  //Libère la mémoire de la séquence
  release_data(seq);
  //Libère la mémoire de la séquence
  release_data(seq2);

  return 0;
}
