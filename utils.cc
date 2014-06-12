#include "utils.h"

void sortKeepAligned(vector<float> &v1, vector<float> &v2) {
    const int n=v1.size();
    /*
    for (int i=0; i<n; i++) cout << v1[i] << " ";
    cout << endl;
    for (int i=0; i<n; i++) cout << v2[i] << " ";
    cout << endl;
    */
    vector<float> w1,w2;
    vector<int> taken;
    for (int i=0; i<n; i++) taken.push_back(0);
    //for (int i=0; i<n; i++) cout << taken[i] << " ";
    //cout << endl;
    
    for (int ii=0; ii<n; ii++) {
        //cout << ii << endl << flush;
        float min=10000000000;
        int ind=-1;
        for (int i=0; i<n; i++) {
            if (v1[i]<=min && !taken[i]) {
                min=v1[i];
                ind=i;
            }
        }
        taken[ind]=1;
        w1.push_back(v1[ind]);
        w2.push_back(v2[ind]);
    }
    v1=w1;
    v2=w2;
}

float defineAverage(vector<float> v) {
    const int n=int(v.size());
    float a=0;
    for (int i=0; i<n; i++) {
        a+=v[i];
    }
    if (n) return a/float(n);
    return 0;
}


//vector<float> defineSmallestIntervalToContain(float interval, vector<float> v) {
float defineSmallestIntervalToContain(float interval, vector<float> v) {
    sort(v.begin(),v.end());
    int nInInterval = int(interval*v.size());
    float minInterval = 100000000000;
    for (int i=0; i<v.size(); i++) {
        //cout << v_m1[i] << " ";
        float interval_size=minInterval;
        if ((i+nInInterval)<v.size()) interval_size = (v[i+nInInterval]-v[i]);
        if (interval_size < minInterval) minInterval = interval_size;
    }
    return minInterval;
}

float invMass(float e1, float eta1, float phi1, float e2, float eta2, float phi2){
    float pt1 = e1/cosh(eta1);
    float pt2 = e2/cosh(eta2);
    float m2 = 2.0*pt1*pt2*(cosh(eta1-eta2)-cos(phi1-phi2));
    return sqrt(m2);
}

float dR(double eta1, double phi1, double eta2, double phi2)
{
  return sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
}

float dRetaphi(double eta1, double phi1, double eta2, double phi2)
{
	if (fabs(phi1-phi2) > pi)
		return sqrt((eta1-eta2)*(eta1-eta2) + pow((2.0*pi - fabs(phi1-phi2)),2.0));
	
	return sqrt((eta1-eta2)*(eta1-eta2) + (phi1-phi2)*(phi1-phi2));
}

float SignalShapeEBfit(double *x, double *par) {
	return SignalShapeEB(x[0]+par[0])*par[1];
}
float SignalShapeEEfit(double *x, double *par) {
	return SignalShapeEE(x[0]+par[0])*par[1];
}
/*
float SignalShapeEBfit(double *x, double *par) {
	return SignalShapeEB(x[0]+par[0])*par[1];
}
*/

float SignalShapeEE(float x) {
	
	const int nBinShape=250;
	float aVec[nBinShape];
	aVec[0] = 9.09091e-05;
	aVec[1] = 9.96016e-05;
	aVec[2] = 8.08219e-05;
	aVec[3] = 0.000119685;
	aVec[4] = 9.56522e-05;
	aVec[5] = 0.000143969;
	aVec[6] = 0.000101639;
	aVec[7] = 0.000115789;
	aVec[8] = 8.31461e-05;
	aVec[9] = 0.000117091;
	aVec[10] = 0.000141224;
	aVec[11] = 5.38745e-05;
	aVec[12] = 0.000158079;
	aVec[13] = 0.000104906;
	aVec[14] = 7.91946e-05;
	aVec[15] = 8e-05;
	aVec[16] = 0.000119703;
	aVec[17] = 9.96516e-05;
	aVec[18] = 8.82591e-05;
	aVec[19] = 0.000146988;
	aVec[20] = 0.000120588;
	aVec[21] = 9.89474e-05;
	aVec[22] = 0.000116788;
	aVec[23] = 0.000151292;
	aVec[24] = 0.000108434;
	aVec[25] = 0.000108651;
	aVec[26] = 8.5259e-05;
	aVec[27] = 0.000146575;
	aVec[28] = 8.66142e-05;
	aVec[29] = 0.000103679;
	aVec[30] = 0.000153307;
	aVec[31] = 7.80328e-05;
	aVec[32] = 8.42105e-05;
	aVec[33] = 0.000114607;
	aVec[34] = 6.90909e-05;
	aVec[35] = 9.22449e-05;
	aVec[36] = 8.92989e-05;
	aVec[37] = 0.000105677;
	aVec[38] = 9.58491e-05;
	aVec[39] = 0.000107383;
	aVec[40] = 8.85714e-05;
	aVec[41] = 0.000119703;
	aVec[42] = 0.000141463;
	aVec[43] = 0.000122267;
	aVec[44] = 6.0241e-05;
	aVec[45] = 0.000120588;
	aVec[46] = 0.000115789;
	aVec[47] = 0.000121168;
	aVec[48] = 8.04428e-05;
	aVec[49] = 9.39759e-05;
	aVec[50] = 3.391e-05;
	aVec[51] = -7.72908e-05;
	aVec[52] = 5.61644e-05;
	aVec[53] = 0.000181102;
	aVec[54] = -3.67893e-05;
	aVec[55] = 6.45914e-05;
	aVec[56] = 2.68852e-05;
	aVec[57] = -1.57895e-05;
	aVec[58] = -2.47191e-05;
	aVec[59] = -5.30909e-05;
	aVec[60] = -1.06122e-05;
	aVec[61] = -3.69004e-06;
	aVec[62] = 5.85153e-05;
	aVec[63] = -3.54717e-05;
	aVec[64] = -9.39597e-06;
	aVec[65] = 5e-05;
	aVec[66] = -3.19703e-05;
	aVec[67] = 3.69338e-05;
	aVec[68] = -0.000169231;
	aVec[69] = -4.09639e-05;
	aVec[70] = -1.61765e-05;
	aVec[71] = -8.21053e-05;
	aVec[72] = 6.86131e-05;
	aVec[73] = -6.12546e-05;
	aVec[74] = -6.50602e-05;
	aVec[75] = 2.56055e-05;
	aVec[76] = -3.90438e-05;
	aVec[77] = -5.89041e-05;
	aVec[78] = 0.000195276;
	aVec[79] = -7.29097e-05;
	aVec[80] = 1.32296e-05;
	aVec[81] = 3.27869e-06;
	aVec[82] = -6.31579e-05;
	aVec[83] = -3.37079e-05;
	aVec[84] = 5.16364e-05;
	aVec[85] = -4.97959e-05;
	aVec[86] = 8.92989e-05;
	aVec[87] = 0.000561572;
	aVec[88] = 0.00171698;
	aVec[89] = 0.0035906;
	aVec[90] = 0.00742143;
	aVec[91] = 0.01249;
	aVec[92] = 0.020922;
	aVec[93] = 0.030336;
	aVec[94] = 0.043588;
	aVec[95] = 0.0593529;
	aVec[96] = 0.0764947;
	aVec[97] = 0.0958978;
	aVec[98] = 0.11567;
	aVec[99] = 0.137087;
	aVec[100] = 0.160643;
	aVec[101] = 0.188973;
	aVec[102] = 0.218062;
	aVec[103] = 0.240682;
	aVec[104] = 0.265773;
	aVec[105] = 0.293327;
	aVec[106] = 0.320832;
	aVec[107] = 0.348689;
	aVec[108] = 0.374411;
	aVec[109] = 0.403575;
	aVec[110] = 0.431147;
	aVec[111] = 0.458644;
	aVec[112] = 0.488217;
	aVec[113] = 0.515254;
	aVec[114] = 0.53902;
	aVec[115] = 0.566819;
	aVec[116] = 0.589752;
	aVec[117] = 0.615587;
	aVec[118] = 0.63808;
	aVec[119] = 0.664137;
	aVec[120] = 0.687959;
	aVec[121] = 0.707379;
	aVec[122] = 0.726826;
	aVec[123] = 0.747043;
	aVec[124] = 0.762952;
	aVec[125] = 0.781342;
	aVec[126] = 0.799515;
	aVec[127] = 0.818559;
	aVec[128] = 0.832301;
	aVec[129] = 0.846658;
	aVec[130] = 0.860791;
	aVec[131] = 0.874009;
	aVec[132] = 0.886242;
	aVec[133] = 0.896982;
	aVec[134] = 0.908713;
	aVec[135] = 0.9189;
	aVec[136] = 0.9283;
	aVec[137] = 0.938453;
	aVec[138] = 0.945995;
	aVec[139] = 0.95186;
	aVec[140] = 0.958876;
	aVec[141] = 0.963283;
	aVec[142] = 0.968989;
	aVec[143] = 0.972846;
	aVec[144] = 0.977824;
	aVec[145] = 0.981341;
	aVec[146] = 0.982743;
	aVec[147] = 0.984082;
	aVec[148] = 0.986609;
	aVec[149] = 0.985925;
	aVec[150] = 0.984707;
	aVec[151] = 0.981384;
	aVec[152] = 0.980941;
	aVec[153] = 0.980132;
	aVec[154] = 0.978176;
	aVec[155] = 0.975987;
	aVec[156] = 0.973396;
	aVec[157] = 0.970395;
	aVec[158] = 0.967422;
	aVec[159] = 0.96346;
	aVec[160] = 0.959622;
	aVec[161] = 0.955396;
	aVec[162] = 0.950846;
	aVec[163] = 0.945995;
	aVec[164] = 0.940835;
	aVec[165] = 0.935373;
	aVec[166] = 0.929544;
	aVec[167] = 0.923836;
	aVec[168] = 0.91787;
	aVec[169] = 0.910812;
	aVec[170] = 0.903924;
	aVec[171] = 0.897669;
	aVec[172] = 0.89127;
	aVec[173] = 0.884613;
	aVec[174] = 0.878051;
	aVec[175] = 0.867734;
	aVec[176] = 0.856756;
	aVec[177] = 0.848477;
	aVec[178] = 0.842189;
	aVec[179] = 0.833787;
	aVec[180] = 0.825903;
	aVec[181] = 0.817931;
	aVec[182] = 0.8102;
	aVec[183] = 0.802771;
	aVec[184] = 0.79388;
	aVec[185] = 0.785348;
	aVec[186] = 0.777295;
	aVec[187] = 0.767854;
	aVec[188] = 0.759619;
	aVec[189] = 0.752358;
	aVec[190] = 0.743394;
	aVec[191] = 0.736674;
	aVec[192] = 0.726982;
	aVec[193] = 0.719321;
	aVec[194] = 0.708937;
	aVec[195] = 0.700003;
	aVec[196] = 0.693202;
	aVec[197] = 0.685961;
	aVec[198] = 0.676889;
	aVec[199] = 0.670513;
	aVec[200] = 0.658788;
	aVec[201] = 0.648343;
	aVec[202] = 0.639381;
	aVec[203] = 0.632761;
	aVec[204] = 0.623614;
	aVec[205] = 0.616164;
	aVec[206] = 0.607884;
	aVec[207] = 0.601437;
	aVec[208] = 0.593944;
	aVec[209] = 0.585817;
	aVec[210] = 0.577248;
	aVec[211] = 0.570328;
	aVec[212] = 0.56147;
	aVec[213] = 0.553509;
	aVec[214] = 0.548765;
	aVec[215] = 0.539964;
	aVec[216] = 0.53343;
	aVec[217] = 0.525787;
	aVec[218] = 0.519868;
	aVec[219] = 0.510817;
	aVec[220] = 0.503274;
	aVec[221] = 0.498478;
	aVec[222] = 0.493058;
	aVec[223] = 0.484176;
	aVec[224] = 0.479896;
	aVec[225] = 0.469998;
	aVec[226] = 0.462893;
	aVec[227] = 0.45576;
	aVec[228] = 0.450375;
	aVec[229] = 0.442683;
	aVec[230] = 0.437887;
	aVec[231] = 0.430736;
	aVec[232] = 0.427216;
	aVec[233] = 0.420775;
	aVec[234] = 0.41485;
	aVec[235] = 0.407769;
	aVec[236] = 0.403546;
	aVec[237] = 0.397028;
	aVec[238] = 0.389534;
	aVec[239] = 0.388569;
	aVec[240] = 0.380296;
	aVec[241] = 0.375312;
	aVec[242] = 0.369457;
	aVec[243] = 0.365884;
	aVec[244] = 0.358889;
	aVec[245] = 0.353344;
	aVec[246] = 0.350364;
	aVec[247] = 0.347136;
	aVec[248] = 0.33872;
	aVec[249] = 0.336706;
	
	int imin=floor(x);
	float y=0;
	if (imin==nBinShape-1) y=aVec[249];
	else {
		y = aVec[imin]+(aVec[imin+1]-aVec[imin])*(x-floor(x)); 
	}
	return y;	
}

float SignalShapeEB(float x) {

	const int nBinShape=250;
	float aVec[nBinShape];
	aVec[0] = 6.94068e-05 ; 
	aVec[1] = -5.03304e-05 ; 
	aVec[2] = -2.13404e-05 ; 
	aVec[3] = 6.017e-05 ; 
	aVec[4] = 2.01697e-05 ; 
	aVec[5] = 0.000114845 ; 
	aVec[6] = 2.13998e-05 ; 
	aVec[7] = 2.74476e-05 ; 
	aVec[8] = 5.2824e-05 ; 
	aVec[9] = 8.754e-05 ; 
	aVec[10] = 2.95346e-06 ; 
	aVec[11] = -7.58699e-05 ; 
	aVec[12] = -2.72224e-05 ; 
	aVec[13] = 3.10997e-06 ; 
	aVec[14] = -3.97771e-05 ; 
	aVec[15] = -1.06916e-05 ; 
	aVec[16] = -0.000113865 ; 
	aVec[17] = 6.05044e-05 ; 
	aVec[18] = -5.81202e-05 ; 
	aVec[19] = -6.58974e-06 ; 
	aVec[20] = 5.37494e-05 ; 
	aVec[21] = -0.000123729 ; 
	aVec[22] = 7.50938e-06 ; 
	aVec[23] = -1.35628e-05 ; 
	aVec[24] = 8.33725e-05 ; 
	aVec[25] = 3.19299e-05 ; 
	aVec[26] = -3.09232e-05 ; 
	aVec[27] = -7.0086e-05 ; 
	aVec[28] = 1.78937e-06 ; 
	aVec[29] = -2.20365e-05 ; 
	aVec[30] = 7.68054e-05 ; 
	aVec[31] = -2.5368e-05 ; 
	aVec[32] = 5.67291e-06 ; 
	aVec[33] = 5.87096e-05 ; 
	aVec[34] = -2.62771e-06 ; 
	aVec[35] = 4.31832e-05 ; 
	aVec[36] = 8.33616e-06 ; 
	aVec[37] = 7.27813e-05 ; 
	aVec[38] = 7.6159e-05 ; 
	aVec[39] = -1.60446e-05 ; 
	aVec[40] = -4.12127e-06 ; 
	aVec[41] = -5.93381e-05 ; 
	aVec[42] = 1.61444e-05 ; 
	aVec[43] = -5.49559e-05 ; 
	aVec[44] = 5.55254e-05 ; 
	aVec[45] = 3.32251e-05 ; 
	aVec[46] = -3.15897e-05 ; 
	aVec[47] = 7.86588e-05 ; 
	aVec[48] = -2.9704e-05 ; 
	aVec[49] = 5.66838e-05 ; 
	aVec[50] = 2.85281e-05 ; 
	aVec[51] = -3.02436e-05 ; 
	aVec[52] = -4.16265e-05 ; 
	aVec[53] = -1.63191e-05 ; 
	aVec[54] = 6.61193e-05 ; 
	aVec[55] = 9.23766e-05 ; 
	aVec[56] = 6.68903e-05 ; 
	aVec[57] = -3.20994e-05 ; 
	aVec[58] = 0.00011082 ; 
	aVec[59] = -4.07997e-05 ; 
	aVec[60] = -8.29046e-06 ; 
	aVec[61] = -7.42197e-05 ; 
	aVec[62] = -1.64386e-05 ; 
	aVec[63] = 1.02508e-05 ; 
	aVec[64] = 7.10995e-06 ; 
	aVec[65] = -5.87486e-05 ; 
	aVec[66] = -0.000101201 ; 
	aVec[67] = 1.62003e-05 ; 
	aVec[68] = -2.53093e-05 ; 
	aVec[69] = 2.65239e-05 ; 
	aVec[70] = -2.68722e-05 ; 
	aVec[71] = -4.02001e-05 ; 
	aVec[72] = 5.0674e-05 ; 
	aVec[73] = -1.75884e-05 ; 
	aVec[74] = 4.7902e-05 ; 
	aVec[75] = -1.01079e-05 ; 
	aVec[76] = 1.08427e-05 ; 
	aVec[77] = -0.000112906 ; 
	aVec[78] = 3.33076e-05 ; 
	aVec[79] = 0.000181201 ; 
	aVec[80] = 0.000426875 ; 
	aVec[81] = 0.00114222 ; 
	aVec[82] = 0.00237804 ; 
	aVec[83] = 0.00541858 ; 
	aVec[84] = 0.0089021 ; 
	aVec[85] = 0.0149157 ; 
	aVec[86] = 0.0231397 ; 
	aVec[87] = 0.0344671 ; 
	aVec[88] = 0.0471013 ; 
	aVec[89] = 0.0625517 ; 
	aVec[90] = 0.0857351 ; 
	aVec[91] = 0.108561 ; 
	aVec[92] = 0.133481 ; 
	aVec[93] = 0.163557 ; 
	aVec[94] = 0.200243 ; 
	aVec[95] = 0.225919 ; 
	aVec[96] = 0.269213 ; 
	aVec[97] = 0.302929 ; 
	aVec[98] = 0.342722 ; 
	aVec[99] = 0.378522 ; 
	aVec[100] = 0.436563 ; 
	aVec[101] = 0.467581 ; 
	aVec[102] = 0.510133 ; 
	aVec[103] = 0.550063 ; 
	aVec[104] = 0.583509 ; 
	aVec[105] = 0.619187 ; 
	aVec[106] = 0.653245 ; 
	aVec[107] = 0.686101 ; 
	aVec[108] = 0.721178 ; 
	aVec[109] = 0.745129 ; 
	aVec[110] = 0.774163 ; 
	aVec[111] = 0.799011 ; 
	aVec[112] = 0.822177 ; 
	aVec[113] = 0.838315 ; 
	aVec[114] = 0.858847 ; 
	aVec[115] = 0.875559 ; 
	aVec[116] = 0.891294 ; 
	aVec[117] = 0.90537 ; 
	aVec[118] = 0.919617 ; 
	aVec[119] = 0.930632 ; 
	aVec[120] = 0.936216 ; 
	aVec[121] = 0.947739 ; 
	aVec[122] = 0.955306 ; 
	aVec[123] = 0.961876 ; 
	aVec[124] = 0.968124 ; 
	aVec[125] = 0.97327 ; 
	aVec[126] = 0.977513 ; 
	aVec[127] = 0.984885 ; 
	aVec[128] = 0.986497 ; 
	aVec[129] = 0.990039 ; 
	aVec[130] = 0.994798 ; 
	aVec[131] = 0.994884 ; 
	aVec[132] = 0.99795 ; 
	aVec[133] = 0.99834 ; 
	aVec[134] = 0.999607 ; 
	aVec[135] = 1 ; 
	aVec[136] = 0.999047 ; 
	aVec[137] = 0.998745 ; 
	aVec[138] = 0.999219 ; 
	aVec[139] = 0.99814 ; 
	aVec[140] = 0.995082 ; 
	aVec[141] = 0.992449 ; 
	aVec[142] = 0.990418 ; 
	aVec[143] = 0.985032 ; 
	aVec[144] = 0.982308 ; 
	aVec[145] = 0.978696 ; 
	aVec[146] = 0.975656 ; 
	aVec[147] = 0.971027 ; 
	aVec[148] = 0.964811 ; 
	aVec[149] = 0.959428 ; 
	aVec[150] = 0.95096 ; 
	aVec[151] = 0.947428 ; 
	aVec[152] = 0.9419 ; 
	aVec[153] = 0.933223 ; 
	aVec[154] = 0.926482 ; 
	aVec[155] = 0.922172 ; 
	aVec[156] = 0.912777 ; 
	aVec[157] = 0.907388 ; 
	aVec[158] = 0.897289 ; 
	aVec[159] = 0.891889 ; 
	aVec[160] = 0.882056 ; 
	aVec[161] = 0.873382 ; 
	aVec[162] = 0.865442 ; 
	aVec[163] = 0.860032 ; 
	aVec[164] = 0.85202 ; 
	aVec[165] = 0.841013 ; 
	aVec[166] = 0.833802 ; 
	aVec[167] = 0.825259 ; 
	aVec[168] = 0.815013 ; 
	aVec[169] = 0.807465 ; 
	aVec[170] = 0.799428 ; 
	aVec[171] = 0.792165 ; 
	aVec[172] = 0.783088 ; 
	aVec[173] = 0.773392 ; 
	aVec[174] = 0.764982 ; 
	aVec[175] = 0.752174 ; 
	aVec[176] = 0.746487 ; 
	aVec[177] = 0.737678 ; 
	aVec[178] = 0.727396 ; 
	aVec[179] = 0.718692 ; 
	aVec[180] = 0.712737 ; 
	aVec[181] = 0.702738 ; 
	aVec[182] = 0.69559 ; 
	aVec[183] = 0.684389 ; 
	aVec[184] = 0.677989 ; 
	aVec[185] = 0.667643 ; 
	aVec[186] = 0.659009 ; 
	aVec[187] = 0.650217 ; 
	aVec[188] = 0.644479 ; 
	aVec[189] = 0.636017 ; 
	aVec[190] = 0.625257 ; 
	aVec[191] = 0.618507 ; 
	aVec[192] = 0.609798 ; 
	aVec[193] = 0.600097 ; 
	aVec[194] = 0.592788 ; 
	aVec[195] = 0.584895 ; 
	aVec[196] = 0.578228 ; 
	aVec[197] = 0.569299 ; 
	aVec[198] = 0.560576 ; 
	aVec[199] = 0.552404 ; 
	aVec[200] = 0.541405 ; 
	aVec[201] = 0.536271 ; 
	aVec[202] = 0.528734 ; 
	aVec[203] = 0.519813 ; 
	aVec[204] = 0.512264 ; 
	aVec[205] = 0.507001 ; 
	aVec[206] = 0.49828 ; 
	aVec[207] = 0.492416 ; 
	aVec[208] = 0.483181 ; 
	aVec[209] = 0.477907 ; 
	aVec[210] = 0.469623 ; 
	aVec[211] = 0.462528 ; 
	aVec[212] = 0.455099 ; 
	aVec[213] = 0.45055 ; 
	aVec[214] = 0.443576 ; 
	aVec[215] = 0.435364 ; 
	aVec[216] = 0.429789 ; 
	aVec[217] = 0.422724 ; 
	aVec[218] = 0.415621 ; 
	aVec[219] = 0.409469 ; 
	aVec[220] = 0.40401 ; 
	aVec[221] = 0.398121 ; 
	aVec[222] = 0.391079 ; 
	aVec[223] = 0.384414 ; 
	aVec[224] = 0.378214 ; 
	aVec[225] = 0.369851 ; 
	aVec[226] = 0.365966 ; 
	aVec[227] = 0.359865 ; 
	aVec[228] = 0.353505 ; 
	aVec[229] = 0.347899 ; 
	aVec[230] = 0.343829 ; 
	aVec[231] = 0.337585 ; 
	aVec[232] = 0.333089 ; 
	aVec[233] = 0.326289 ; 
	aVec[234] = 0.322249 ; 
	aVec[235] = 0.316079 ; 
	aVec[236] = 0.31061 ; 
	aVec[237] = 0.305426 ; 
	aVec[238] = 0.301885 ; 
	aVec[239] = 0.296753 ; 
	aVec[240] = 0.290931 ; 
	aVec[241] = 0.286877 ; 
	aVec[242] = 0.281831 ; 
	aVec[243] = 0.276633 ; 
	aVec[244] = 0.272283 ; 
	aVec[245] = 0.268069 ; 
	aVec[246] = 0.26399 ; 
	aVec[247] = 0.258457 ; 
	aVec[248] = 0.253549 ; 
	aVec[249] = 0.249493 ; 
	
	int imin=floor(x);
	float y=0;
	if (imin==nBinShape-1) y=aVec[249];
	else {
		y = aVec[imin]+(aVec[imin+1]-aVec[imin])*(x-floor(x)); 
	}
	return y;
}