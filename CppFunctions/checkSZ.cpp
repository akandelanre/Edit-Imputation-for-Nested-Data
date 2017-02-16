#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;
  
#define HEAD 1
#define SPOUSE 2
#define BIOLOGICALCHILD 3
#define ADOPTEDCHILD 4
#define STEPCHILD 5
#define SIBLING 6
#define PARENT 7
#define GRANDCHILD 8
#define PARENTINLAW 9
#define CHILDINLAW 10

#define GENDER 0
#define AGE 3
#define RELATE 4

inline int GetHead(NumericMatrix hh_to_check, int h) {
  for(int kk = 0; kk < h; kk++){
    if (hh_to_check(kk,RELATE) == HEAD) {
      return kk;
    }
  }
  return -1;
}

inline bool IsHead(double relate, double age) {
  return (relate == HEAD && age >=16);
}

inline bool MoreThanOneHead(NumericMatrix hh_to_check, int h) {
  int nhead = 0;
  for(int kk = 0; kk < h; kk++){
    if (hh_to_check(kk,RELATE) == HEAD) {
      nhead++;
    }
  }
  return (nhead >1);
}

inline int GetValidSpouse(NumericMatrix hh_to_check, int h) {
  int nspouse = 0;
  int spouse = -1;
  for(int kk = 0; kk < h; kk++){
    if (hh_to_check(kk,RELATE) == SPOUSE) {
      nspouse++;
      spouse = kk;
    }
  }
  if (nspouse > 1) {return -2;} //too many spouse
  if (nspouse == 0) { return -1;} //no spouse
  if (hh_to_check(spouse,AGE)<16) {return -2;} //spouse is under-age
  return spouse;
}

inline bool IsValidCouple(NumericMatrix hh_to_check, int spouse, int head) {
  if (spouse == -2) { //bad spouse or too many spouses
    return false;
  } else { //valid spouse or no spouse
    if (spouse >= 0) {//the only spouse, so check sex, and age difference
      if (hh_to_check(head,GENDER) == hh_to_check(spouse,GENDER)) {return false;}
      if (std::abs(hh_to_check(head,AGE) - hh_to_check(spouse,AGE)) > 49) {return false;}
    }
  }
  return true;
}

//return -1 if no biological child
//return the record index of the oldest biological child otherwise
inline int GetOldestBiologicalChild(NumericMatrix hh_to_check, int h) {
  double age = -1;
  int child = -1;  //no biological childen
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)==BIOLOGICALCHILD) {
      if (hh_to_check(i,AGE) > age) {
        age = hh_to_check(i,AGE);
        child = i;
      }
    }
  }
  return child;
}

inline bool IsValidBiologicalChild(NumericMatrix hh_to_check, int child, int head) {
  if (child>=0) {//get a biological child, check age difference
    if (hh_to_check(head,AGE) - hh_to_check(child,AGE) < 7) {return false;}
  }
  return true;
}


//return -1 if no adopted child
//return the record index of the oldest adopted child otherwise
inline int GetOldestAdoptedChild(NumericMatrix hh_to_check, int h) {
  double age = -1;
  int child = -1;  //no adopted childen
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)==ADOPTEDCHILD) {
      if (hh_to_check(i,AGE) > age) {
        age = hh_to_check(i,AGE);
        child = i;
      }
    }
  }
  return child;
}

inline bool IsValidAdoptedChild(NumericMatrix hh_to_check, int child, int head) {
  if (child>=0) {//get an adopted child, check age difference
    if (hh_to_check(head,AGE) - hh_to_check(child,AGE) <11) {return false;}
  }
  return true;
}


//return -1 if no step child
//return the record index of the oldest step child otherwise
inline int GetOldestStepChild(NumericMatrix hh_to_check, int h) {
  double age = -1;
  int child = -1;  //no step childen
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)==STEPCHILD) {
      if (hh_to_check(i,AGE) > age) {
        age = hh_to_check(i,AGE);
        child = i;
      }
    }
  }
  return child;
}

inline bool IsValidStepChild(NumericMatrix hh_to_check, int child, int head) {
  if (child>=0) {//get a step child, check age difference
    if (hh_to_check(head,AGE) - hh_to_check(child,AGE) <9) {return false;}
  }
  return true;
}


//return -1 if no parent
//return the record index of the youngest parent otherwise
inline int GetYoungestParent(NumericMatrix hh_to_check, int h) {
  double age = 1000;
  int parent = -1;  //no parent
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)==PARENT) {
      if (hh_to_check(i,AGE) < age) {
        age = hh_to_check(i,AGE);
        parent = i;
      }
    }
  }
  return parent;
}

inline bool IsValidParent(NumericMatrix hh_to_check, int parent, int head) {
  if (parent >= 0) {//get a child, check age difference
    if (hh_to_check(parent,AGE) - hh_to_check(head,AGE) < 4) {return false;}
  }
  return true;
}

//return -1 if no parentinlaw
//return the record index of the youngest parentinlaw otherwise
inline int GetYoungestParentInLaw(NumericMatrix hh_to_check, int h) {
  double age = 1000;
  int parent = -1;  //no parent
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)==PARENTINLAW) {
      if (hh_to_check(i,AGE) < age) {
        age = hh_to_check(i,AGE);
        parent = i;
      }
    }
  }
  return parent;
}

inline bool IsValidParentInLaw(NumericMatrix hh_to_check, int parentinlaw, int head) {
  if (parentinlaw >= 0) {//get a child, check age difference
    if (hh_to_check(parentinlaw,AGE) - hh_to_check(head,AGE) < 4) {return false;}
  }
  return true;
}


inline bool IsValidSibling(NumericMatrix hh_to_check, int h, int head) {
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE) == SIBLING) {
      if (std::abs(hh_to_check(i,AGE) - hh_to_check(head,AGE)) > 37) {return false;}
    }
  }
  return true;
}

inline bool IsValidGrandChild(NumericMatrix hh_to_check, int h, int spouse, int head) {
  for (int i = 0; i < h; i++) {
    if (hh_to_check(i,RELATE)== GRANDCHILD) {
      if (hh_to_check(head,AGE) < 31) {return false;} //too young to be grand parent for the HEAD
      if (spouse >= 0) { //make sure the spouse(if any) is not too young
        if (hh_to_check(spouse,AGE) < 17) {return false;}
      }
      if (hh_to_check(head,AGE) - hh_to_check(i,AGE) < 26) {return false;}
    }
  }
  return true;
}

inline int isValid(NumericMatrix hh_to_check, int h) {
  
  int head = GetHead(hh_to_check,h);
  if (head < 0) {return 0;}
  
  if (!IsHead(hh_to_check(head,RELATE), hh_to_check(head,AGE))) {return 0;}
  if (MoreThanOneHead(hh_to_check,h)) {return 0;}
  
  int spouse = GetValidSpouse(hh_to_check,h);
  if (!IsValidCouple(hh_to_check,spouse, head)) {return 0;}
  
  int oldestBiologicalChild = GetOldestBiologicalChild(hh_to_check,h);
  if (!IsValidBiologicalChild(hh_to_check,oldestBiologicalChild,head)) {return 0;}
  
  int oldestAdoptedChild = GetOldestAdoptedChild(hh_to_check,h);
  if (!IsValidAdoptedChild(hh_to_check,oldestAdoptedChild,head)) {return 0;}
  
  int oldestStepChild = GetOldestStepChild(hh_to_check,h);
  if (!IsValidStepChild(hh_to_check,oldestStepChild,head)) {return 0;}
  
  int youngestParent = GetYoungestParent(hh_to_check,h);
  if (!IsValidParent(hh_to_check,youngestParent,head)) {return 0;}
  
  int youngestParentInLaw = GetYoungestParentInLaw(hh_to_check,h);
  if (!IsValidParentInLaw(hh_to_check,youngestParentInLaw,head)) {return 0;}
  
  if (!IsValidSibling(hh_to_check,h,head)) {return 0;}
  
  if (!IsValidGrandChild(hh_to_check,h,spouse,head)) {return 0;}
  
  return 1;
}

// [[Rcpp::export]]
NumericVector checkSZ(NumericMatrix Data_to_check, int h){
  int n0 = Data_to_check.nrow();
  int p = Data_to_check.ncol()/h;
  
  NumericVector Data_checked(n0);
  
  for (int i = 0; i < n0; i++) {
    NumericMatrix hh_to_check(h, p);
    for(int j = 0; j < p; j++){
      for(int k = 0; k < h; k++){
        hh_to_check(k,j) = Data_to_check(i,j+(k*p));
      }
    }
    Data_checked[i] = isValid(hh_to_check,h);
  }
  return Data_checked;
}




