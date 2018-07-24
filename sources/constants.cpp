#include "constants.h"

using namespace std;


const vector<int> Const_arrays::belt( [](){
                                          vector<int> belt(500);
                                          belt[0] = 1;
                                          for(unsigned i=1; i< belt.size(); i++) belt[i] = -belt[i-1];
                                          return belt;
                                      }()
                                      );

const vector<double> Const_arrays::fact( [](){
                                             vector<double> fact(160);
                                             fact[0] = 1.0;
                                             fact[1] = 1.0;
                                             for(int i=2; i< fact.size(); i++) fact[i] = fact[i-1] * i;
                                             return fact;
                                         }()
                                         );

const vector<double> Const_arrays::dfact( [](){
                                              vector<double> dfact(200);
                                              dfact[0] = 1.0;
                                              dfact[1] = 1.0;
                                              for(unsigned i = 1; i < dfact.size(); i++) dfact[i+1] = dfact[i] * ( 2*i + 1 );
                                              return dfact;
                                          }()
                                          );

const vector<vector<double>> Const_arrays::binom( [](){
                                                      vector<vector<double>> binom(100);
                                                      for(unsigned i = 0; i < binom.size(); i++) {
                                                          binom[i].resize(binom.size());
                                                          for(unsigned j = 0; j <= i; j++) {
                                                              binom[i][j] = Const_arrays::fact[i];
                                                              binom[i][j] /= Const_arrays::fact[j] * Const_arrays::fact[i-j];
                                                              binom[j][i] = binom[i][j];
                                                          }
                                                      }
                                                      return binom;
                                                  }()
                                                  );

const vector<vector<double>> Const_arrays::omega( []() {
                                                      vector<vector<double>> omega(33);
                                                      for(unsigned i = 0; i <= 32; i++) {
                                                          omega[i].resize(33);
                                                          for(unsigned j = 0; j <= i; j++)
                                                          omega[i][j] = std::sqrt( ( 2.0*i + 1 ) * Const_arrays::fact[i-j] / ( 2.0 * Const_arrays::fact[i+j] ) );
                                                      }
                                                      return omega;
                                                  }()
                                                  );
