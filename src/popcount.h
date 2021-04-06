//Compte le nombre de bits a 1 dans l'octet
int popcount(char c){
  int pop=0;
  //On fait 8 tours pour parcourir tout l'octet, bit par bit
  for(int i=0;i<8;i++){
    //On ajoute 0 ou 1 selon le LSB (Least Significant Digit), bit de pois faible
    pop+= c & 1;
    //On shift a droite
    /*
    c : 0100 0101 
    c = c >> 1
    c : 0010 0010
    */
    c = c >> 1;
  }
  return pop;
}
