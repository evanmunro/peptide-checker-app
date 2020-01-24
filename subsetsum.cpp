#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
using namespace std;
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


double cap_mass = 41.0265;
double WATER_MASS = 18.0106;
double tolerance = 1.0; 
CharacterVector matched_peptides(1); 
NumericVector matched_masses(1); 


map<char, double> AMINO_MASS = {
  { 'A', 89.0477 },
  { 'R', 174.1117 },
  { 'N', 132.0535},
  {'D', 133.0375},
  { 'C', 121.0197},
  { 'E', 147.0532 },
  { 'Q', 146.0691 },
  {'G', 75.032},
  { 'H', 155.0695 },
  { 'I', 131.0946 },
  { 'L', 131.0946},
  {'K', 146.1055},
  { 'M', 149.051 },
  { 'F', 165.079 },
  { 'P', 115.0633},
  {'S', 105.0426},
  { 'T', 119.0582},
  { 'W', 204.0899 },
  { 'Y', 181.0739 },
  {'V', 117.079}
};

double sum_vector(vector<double>& masses) {
  double sum  = 0.0;
  for (int i = 0 ; i < masses.size(); i++){
    sum+= masses[i];
  }
  return sum;
}

//recursive function to search for side products of a peptide synthesis of
//a certain mass
void findCandidatePeptides(vector<char>& orig_pep,double target_mass,vector<char>& component_aminos,
                           vector<double>& component_masses, int end_ignore, int n) {
  double candidate_mass = sum_vector(component_masses);
  candidate_mass = candidate_mass - WATER_MASS*(component_masses.size()-1) + cap_mass;
  
  if(abs(candidate_mass-target_mass)<tolerance){
    matched_masses.push_back(candidate_mass); 
    vector<char> side_product(orig_pep.begin(),orig_pep.end());
    int j = orig_pep.size()-1;
    int i = component_aminos.size()-1;
    while(j>=0 && i >=0 ) {
      if(component_aminos[i]== orig_pep[j]){
        j= j-1;
        i=i-1;
      }
      else{
        side_product[j] = '*';
        j= j-1;
      }
    }
    while(j>=0) {
      side_product[j] = '*';
      j= j-1;
    }
    string s(side_product.begin(),side_product.end());
    matched_peptides.push_back(s); 
  }
  
  if(candidate_mass < target_mass) {
    return;
  }
  for (int i = n ; i < component_masses.size()-end_ignore; i++) {
    vector<double> drop_one_mass(component_masses.size());
    copy(component_masses.begin(), component_masses.end(),drop_one_mass.begin());
    vector<char> drop_one_acids(component_aminos.size());
    copy(component_aminos.begin(), component_aminos.end(),drop_one_acids.begin());
    drop_one_mass.erase(drop_one_mass.begin()+i);
    drop_one_acids.erase(drop_one_acids.begin()+i);
    findCandidatePeptides(orig_pep,target_mass,drop_one_acids,drop_one_mass,end_ignore,i);
  }
  
}

//[[Rcpp::export]]
List searchPeptides(string peptide, NumericVector adjustments,double mass, double tol, int ignore,double cap) { 
  cap_mass = cap; 
  matched_peptides = CharacterVector(1); 
  matched_masses = NumericVector(1); 
  tolerance = tol; 
  vector<char> pept(peptide.begin(),peptide.end());
  vector<double> arr(pept.size());
  for (int i = 0; i < (pept.size()); i++) {
    arr[i] = AMINO_MASS[pept[i]]+adjustments[i];
  }
  findCandidatePeptides(pept,mass,pept,arr,ignore,0); 
  matched_peptides.erase(0); 
  matched_masses.erase(0); 
  List result = List::create(Named("Peptide")=matched_peptides,Named("Mass")=matched_masses);
  return(result); 
}
// [[Rcpp::export]]
List listTruncated(string peptide,NumericVector adjustments, double cap) { 
  cap_mass = cap; 
  vector<char> pept(peptide.begin(),peptide.end());
  vector<double> arr(pept.size());
  for (int i = 0; i < (pept.size()); i++) {
    arr[i] = AMINO_MASS[pept[i]]+adjustments[i];
  }
  int n_acids = pept.size()-1; 
  CharacterVector truncated(n_acids); 
  NumericVector mass(n_acids); 
  List result = List::create(Named("peptide")=truncated,Named("mass")=mass); 
  for (int k=1; k <9; k++) {
    NumericVector v(n_acids); 
    result.push_back(v,"k="+to_string(k)); 
  }
  
  int row=0; 
  double running_sum = arr[pept.size()-1];
  string running_pept= peptide.substr(pept.size()-1,1);
  for (int i =pept.size()-2; i >=0 ; i=i-1) {
    running_sum += arr[i]-WATER_MASS;
    running_pept = peptide[i] + running_pept;
    truncated(row)= running_pept; 
    mass(row) = running_sum + cap_mass; 

    for (int k =1; k < 9; k++) {
      NumericVector v = result[k+1]; 
      v(row) = (running_sum+cap_mass+k)/k; 
    }
    row++; 
  }
  return(result); 
}

