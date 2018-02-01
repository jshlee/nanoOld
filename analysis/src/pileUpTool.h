#include <vector>
#include <numeric>

class pileUpTool {

  std::vector<float> m_weights;
  std::vector<float> m_weights_up;
  std::vector<float> m_weights_dn;

  static const std::vector<float> Moriond17MC;
  static const std::vector<float> Moriond17RD;
  static const std::vector<float> Moriond17RD_up;
  static const std::vector<float> Moriond17RD_dn;

 public:  
  pileUpTool(){
    m_weights.clear(); m_weights_up.clear(); m_weights_dn.clear();

    std::vector<float> pileupMC = Moriond17MC;
    std::vector<float> pileupRD = Moriond17RD;
    std::vector<float> pileupUp = Moriond17RD_up;
    std::vector<float> pileupDn = Moriond17RD_dn;
    
    const double sumWMC = std::accumulate(pileupMC.begin(), pileupMC.end(), 0.);
    const double sumWRD = std::accumulate(pileupRD.begin(), pileupRD.end(), 0.);
    const double sumWUp = std::accumulate(pileupUp.begin(), pileupUp.end(), 0.);
    const double sumWDn = std::accumulate(pileupDn.begin(), pileupDn.end(), 0.);

    std::vector<float> pileupMCTmp;
    std::vector<float> pileupRDTmp;
    std::vector<float> pileupUpTmp, pileupDnTmp;

    for ( int i=0, n=std::min(pileupMC.size(), pileupRD.size()); i<n; ++i ){
      pileupMCTmp.push_back(pileupMC[i]/sumWMC);
      pileupRDTmp.push_back(pileupRD[i]/sumWRD);
      pileupUpTmp.push_back(pileupUp[i]/sumWUp);
      pileupDnTmp.push_back(pileupDn[i]/sumWDn);
    }
    
    for ( int i=0, n=std::min(pileupMC.size(), pileupRD.size()); i<n; ++i ){
      m_weights.push_back(pileupRDTmp[i]/pileupMCTmp[i]);
      m_weights_up.push_back(pileupUpTmp[i]/pileupMCTmp[i]);
      m_weights_dn.push_back(pileupDnTmp[i]/pileupMCTmp[i]);
      /* std::cout << " " << i */
      /* 		<< " weight = "<< pileupRDTmp[i]/pileupMCTmp[i] */
      /* 		<< std::endl; */
    }
  };
  
  float getWeight(int nTrueInt, int sys = 0){
    if (sys == 0) return m_weights[nTrueInt];
    if (sys == 1) return m_weights_up[nTrueInt];
    if (sys == -1) return m_weights_dn[nTrueInt];
    return 1.0;
  };
  
};

const std::vector<float> pileUpTool::Moriond17MC = {
  1.78653e-05,2.56602e-05,5.27857e-05,8.88954e-05,0.000109362,
  0.000140973,0.000240998,0.00071209 ,0.00130121 ,0.00245255 ,
  0.00502589 ,0.00919534 ,0.0146697  ,0.0204126  ,0.0267586  ,
  0.0337697  ,0.0401478  ,0.0450159  ,0.0490577  ,0.0524855  ,
  0.0548159  ,0.0559937  ,0.0554468  ,0.0537687  ,0.0512055  ,
  0.0476713  ,0.0435312  ,0.0393107  ,0.0349812  ,0.0307413  ,
  0.0272425  ,0.0237115  ,0.0208329  ,0.0182459  ,0.0160712  ,
  0.0142498  ,0.012804   ,0.011571   ,0.010547   ,0.00959489 ,
  0.00891718 ,0.00829292 ,0.0076195  ,0.0069806  ,0.0062025  ,
  0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,
  0.00212473 ,0.00160208 ,0.00117884 ,0.000859662,0.000569085,
  0.000365431,0.000243565,0.00015688 ,9.88128e-05,6.53783e-05,
  3.73924e-05,2.61382e-05,2.0307e-05 ,1.73032e-05,1.435e-05  ,
  1.36486e-05,1.35555e-05,1.37491e-05,1.34255e-05,1.33987e-05,
  1.34061e-05,1.34211e-05,1.34177e-05,1.32959e-05,1.33287e-05
};

const std::vector<float> pileUpTool::Moriond17RD = {
  2.387970e+05,8.375429e+05,2.308427e+06,3.124754e+06,4.476191e+06,
  5.995911e+06,7.000896e+06,1.289165e+07,3.526173e+07,7.870123e+07,
  1.769458e+08,3.600895e+08,6.027665e+08,8.765194e+08,1.174474e+09,
  1.489059e+09,1.759352e+09,1.943926e+09,2.049172e+09,2.101582e+09,
  2.132787e+09,2.149099e+09,2.128986e+09,2.062649e+09,1.962884e+09,
  1.841872e+09,1.704136e+09,1.554523e+09,1.399489e+09,1.243533e+09,
  1.088821e+09,9.373048e+08,7.920441e+08,6.567177e+08,5.344668e+08,
  4.271268e+08,3.351056e+08,2.577246e+08,1.937514e+08,1.418309e+08,
  1.006714e+08,6.901386e+07,4.554008e+07,2.884748e+07,1.750632e+07,
  1.016264e+07,5.637781e+06,2.987282e+06,1.512002e+06,7.318454e+05,
  3.398220e+05,1.525454e+05,6.740482e+04,3.048969e+04,1.515211e+04,
  8.975911e+03,6.496155e+03,5.434805e+03,4.889958e+03,4.521716e+03,
  4.208464e+03,3.909763e+03,3.614274e+03,3.320722e+03,3.031096e+03,
  2.748237e+03,2.474977e+03,2.213817e+03,1.966815e+03,1.735546e+03,
  1.521109e+03,1.324149e+03,1.144898e+03,9.832202e+02,8.386676e+02,
};
const std::vector<float> pileUpTool::Moriond17RD_up = {
  2.326834e+05,6.594790e+05,2.183879e+06,2.745838e+06,4.071578e+06,
  5.399946e+06,6.385996e+06,9.041865e+06,2.378316e+07,5.400307e+07,
  1.160154e+08,2.460525e+08,4.433853e+08,6.801119e+08,9.378723e+08,
  1.218825e+09,1.501706e+09,1.730477e+09,1.882658e+09,1.968697e+09,
  2.012675e+09,2.040380e+09,2.054745e+09,2.036917e+09,1.978318e+09,
  1.889436e+09,1.781249e+09,1.657942e+09,1.523136e+09,1.382240e+09,
  1.239666e+09,1.097667e+09,9.577008e+08,8.219249e+08,6.933070e+08,
  5.748163e+08,4.686303e+08,3.757644e+08,2.961701e+08,2.290718e+08,
  1.733771e+08,1.279604e+08,9.175390e+07,6.370411e+07,4.270206e+07,
  2.757163e+07,1.711756e+07,1.020509e+07,5.837087e+06,3.201566e+06,
  1.683875e+06,8.498863e+05,4.125384e+05,1.935936e+05,8.888115e+04,
  4.097021e+04,1.993457e+04,1.101254e+04,7.298604e+03,5.722913e+03,
  4.985459e+03,4.560494e+03,4.244754e+03,3.963531e+03,3.691485e+03,
  3.421626e+03,3.153557e+03,2.889088e+03,2.630613e+03,2.380501e+03,
  2.140853e+03,1.913419e+03,1.699565e+03,1.500273e+03,1.316159e+03,
};
const std::vector<float> pileUpTool::Moriond17RD_dn = {
  2.474113e+05,1.069278e+06,2.428277e+06,3.566992e+06,4.991832e+06,
  6.593141e+06,8.097305e+06,1.996320e+07,5.191520e+07,1.197807e+08,
  2.727752e+08,5.131701e+08,8.023341e+08,1.118763e+09,1.462995e+09,
  1.780219e+09,2.005191e+09,2.135073e+09,2.198487e+09,2.233964e+09,
  2.252561e+09,2.229714e+09,2.154206e+09,2.041681e+09,1.905607e+09,
  1.751084e+09,1.584596e+09,1.413602e+09,1.242528e+09,1.073635e+09,
  9.097710e+08,7.551027e+08,6.138886e+08,4.891245e+08,3.819822e+08,
  2.920568e+08,2.180117e+08,1.582197e+08,1.111016e+08,7.513638e+07,
  4.874432e+07,3.023911e+07,1.789611e+07,1.008679e+07,5.408384e+06,
  2.757256e+06,1.336896e+06,6.175226e+05,2.730083e+05,1.168838e+05,
  4.983506e+04,2.245556e+04,1.173965e+04,7.637138e+03,6.017980e+03,
  5.280664e+03,4.837353e+03,4.483816e+03,4.153851e+03,3.828857e+03,
  3.506027e+03,3.187482e+03,2.876619e+03,2.576851e+03,2.291178e+03,
  2.022034e+03,1.771241e+03,1.540024e+03,1.329043e+03,1.138448e+03,
  9.679458e+02,8.168710e+02,6.842602e+02,5.689249e+02,4.695210e+02,
};
