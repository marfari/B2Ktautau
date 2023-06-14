
using namespace std;

ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag);
ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag);

// Global variables
int dimM = 24; // number of known parameters
int dimX = 23; // number of unkown parameters

TVectorD m( dimM );
TMatrixD V( dimM, dimM );

TVectorD x_initial_estimate( TVectorD m );
TVectorD h( TVectorD x );
TMatrixD dh_dx( TVectorD x );

TVectorD first_chi2_derivative( TVectorD x );
TMatrixD second_chi2_derivative( TVectorD x );

double chisquare( const double* x_values );

// void minimise( TVectorD x );

double mtau = 1776.86; // MeV (PDG 2023)

void decay_fit()
{
    TFile* f = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2018/DVntuple_MC_2018_MagUp_visible.root");
    TTree* t = (TTree*)f->Get("ntuple/DecayTree");

    UInt_t num_entries = t->GetEntries(); 
    Double_t m_vars[dimM];
    Double_t V_vars[dimM][dimM];

    for(int i = 1; i <= dimM; i++)
    {
        t->SetBranchAddress(Form("df_m_%i",i), &m_vars[i-1] );
        for(int j = 1; j <= dimM; j++)
        {
            t->SetBranchAddress(Form("df_V_%i%i",i,j), &V_vars[i-1][j-1]);
        }
    }

    TFile* fout = new TFile("/panfs/felician/B2Ktautau/ROOT_Sim/2018/decay_fit_tuple_2018.root", "RECREATE");
    TTree* tout = new TTree("DecayTree", "DecayTree");

    Double_t X[dimX];
    Double_t XERR[dimX];
    double chi2, MB;
    int status;

    for(int i = 0; i < dimX; i++)
    {
      tout->Branch(Form("df_x_%i",i), &X[i]);
      tout->Branch(Form("df_x_error_%i",i), &XERR[i]);
    }
    tout->Branch("df_chi2", &chi2);
    tout->Branch("df_Mb", &MB);
    tout->Branch("df_status", &status);

    // Loop over events
    for(int evt = 0; evt < 1; evt++)
    {
      t->GetEntry(evt);
      
      // 1)  Measured parameters, m, and their covariance matrix V
      // m = {PV, DV1, P3pi1, DV2, P3pi1, RP, PK}
      for(int i = 0; i < dimM; i++)
      {
          m(i) = m_vars[i];
          for(int j = 0; j < dimM; j++)
          {
            V(i,j) = V_vars[i][j];
          }
      }
      V(19,19) = pow(10,-5); // z-component of the reference point is fixed; its error is zero; i'm setting it to a small non-zero number

      // 2) Get an initial estimate for x = {BV, PB, ptau1, pnu1, ptau2, pnu2, L1, L2, L, LK}
      TVectorD x0 = x_initial_estimate(m);

      // minimise(x0);
  
      ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2");    
      ROOT::Math::Functor func(&chisquare, dimX);

      min->SetMaxIterations(10000000);  
      min->SetPrintLevel(2);
      min->SetPrecision(pow(10,-20));
      
      min->SetFunction(func);

      double x0_values[dimX];
      double x0_err[dimX];
      for(int i = 0; i < dimX; i++)
      {
        x0_values[i] = x0(i);
        x0_err[i] = 0.1*abs(x0_values[i]);

        min->SetVariable(i, Form("x_%i",i), x0_values[i], x0_err[i]);
        // min->SetVariableLimits(i, x_min[i], x_max[i]);
      }

      min->Minimize();

      min->PrintResults();

      chi2 = min->MinValue();
      status = min->Status();

      const double *xMin = min->X();
      for(int i = 0; i < dimX; i++)
      {
        X[i] = xMin[i];
      }

      ROOT::Math::XYZVector pB_min( X[3], X[4], X[5] );
      double EB_min = X[6];
      if( pow(EB_min,2) > pB_min.Mag2() )
      {
        MB = sqrt( pow(EB_min,2) - pB_min.Mag2() );
      }
      else
      {
        MB = -sqrt( pB_min.Mag2() - pow(EB_min,2) );
      }

      cout << "initial chi2 = " << chisquare(x0_values) << endl;
      cout << "final chi2 = " << min->MinValue() << endl;
      cout << "MB = " << MB << endl;
      cout << "Status = " << status << endl;

      tout->Fill();
    }
    fout->cd();
    tout->Write();
    fout->Close();
}

// void minimise( TVectorD x )
// {
//   TVectorD first(dimX);
//   TMatrixD second(dimX,dimX);
//   TMatrixD second_inverse(dimX,dimX);

//   TVectorD dx(dimX);
//   TVectorD x_new(dimX);
//   int max_iterations = 1000;

//   TVectorD x0 = x;

//   for(int i = 0; i < max_iterations; i++)
//   {

//     first = first_chi2_derivative(x);
//     second = second_chi2_derivative(x);
//     second_inverse = second;
//     second_inverse.Invert();

//     dx = second_inverse*first;

//     x_new = x - dx;

//     cout << chisquare(x_new) << endl;

//     if( abs(chisquare(x) - chisquare(x_new)) < 0.001 )
//     {
//       break;
//     }
//     else
//     {
//       x = x_new;
//     }

//     if(i == max_iterations-1)
//     {
//       cout << "Maximum number of iterations has been reached. No convergence." << endl;
//     }
//   }

//   cout << "Initial chi2 = " << chisquare(x0) << endl;
//   cout << "Final chi2 = " << chisquare(x_new) << endl;
// }

double chisquare( const double* x_values )
{
    TVectorD x(dimX);
    for(int i = 0; i < dimX; i++)
    {
      x(i) = x_values[i];
    }

    TVectorD hx = h( x );

    TMatrixD W = V;
    W.Invert();

    //TVectorD c = W*(m - hx);

    double chi2 = (m-hx)*(W*(m - hx));
    return chi2;
}

TMatrixD second_chi2_derivative( TVectorD x )
{
  // here we assume that h(x) varies slowly with x such that its second derivative is ~ 0
  TMatrixD H = dh_dx(x);

  TMatrixD H_transpose = H;
  H_transpose.T();

  TMatrixD V_inverse = V;
  V_inverse.Invert();

  TMatrixD second_derivative = H_transpose*V_inverse*H;
  second_derivative *= 2;

  return second_derivative;
}


TVectorD first_chi2_derivative( TVectorD x )
{
  TVectorD hx = h(x);
  TMatrixD H = dh_dx(x);

  TMatrixD H_transpose = H;
  H_transpose.T();

  TMatrixD V_inverse = V;
  V_inverse.Invert();

  TVectorD first_derivative = H_transpose*V_inverse*(m - hx);
  first_derivative *= -2;

  return first_derivative;

}

TMatrixD dh_dx( TVectorD x ) 
{
  // H = d h(x)/dx

  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  double EB = x(6);
  ROOT::Math::XYZVector ptau1( x(7), x(8), x(9) );
  ROOT::Math::XYZVector pnu1( x(10), x(11), x(12) );
  ROOT::Math::XYZVector ptau2( x(13), x(14), x(15) );
  ROOT::Math::XYZVector pnu2( x(16), x(17), x(18) );
  double L1 = x(19);
  double L2 = x(20);
  double L = x(21);
  double LK = x(22);

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  double Enu1 = sqrt( pnu1.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );

  double h1[] = {	1,	0,	0,	-L,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-pB.x(),	0 };
  double h2[] = {	0,	1,	0,	0,	-L,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-pB.y(),	0 };
  double h3[] = { 0,	0,	1,	0,	0,	-L,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	-pB.z(),	0 };
  double h4[] = { 1,	0,	0,	0,	0,	0,	0,	L1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	ptau1.x(),	0,	0,	0 };
  double h5[] = { 0,	1,	0,	0,	0,	0,	0,	0,	L1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	ptau1.y(),	0,	0,	0 };
  double h6[] = { 0,	0,	1,	0,	0,	0,	0,	0,	0,	L1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	ptau1.z(),	0,	0,	0 };
  double h7[] = { 0,	0,  0,	0,	0,  0,  0,	1,  0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h8[] = { 0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h9[] = { 0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h10[] = { 0,	0,	0,	0,	0,	0,	0,	ptau1.x()/Etau1, 	ptau1.y()/Etau1,	ptau1.z()/Etau1,	-pnu1.x()/Enu1,	-pnu1.y()/Enu1,	-pnu1.z()/Enu1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h11[] = { 1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	L2,	0,	0,	0,	0,	0,	0,	ptau2.x(),	0,	0 };
  double h12[] = { 0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	L2,	0,	0,	0,	0,	0,	ptau2.y(),	0,	0 };
  double h13[] = { 0,	0,	1,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	L2,	0,	0,	0,	0,	ptau2.z(),	0,	0 };
  double h14[] = { 0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	-1,	0,	0,	0,	0,	0,	0 };
  double h15[] = { 0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	-1,	0,	0,	0,	0,	0 };
  double h16[] = { 0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	1,	0,	0,	-1,	0,	0,	0,	0 };
  double h17[] = { 0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	ptau2.x()/Etau2,	ptau2.y()/Etau2,	ptau2.z()/Etau2,	-pnu2.x()/Enu2,	-pnu2.y()/Enu2,	-pnu2.z()/Enu2,	0,	0,	0,	0 };
  double h18[] = { 1,	0,	0,	LK,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	0,	0,	0,	pB.x() - ptau1.x() - ptau2.x() };
  double h19[] = { 0,	1,	0,	0,	LK,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	0,	0,	pB.y() - ptau1.y() - ptau2.y() };
  double h20[] = { 0,	0,	1,	0,	0,	LK,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	-LK,	0,	0,	0,	0,	0,	0,	pB.z() - ptau1.z() - ptau2.z() };
  double h21[] = { 0,	0,	0,	1,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h22[] = { 0,	0,	0,	0,	1,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0,	0 };
  double h23[] = { 0,	0,	0,	0,	0,	1,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	-1,	0,	0,	0,	0,	0,	0,	0 };
  double h24[] = { 0,	0,	0,	0,	0,	0,	1,	-ptau1.x()/Etau1,	-ptau1.y()/Etau1,	-ptau1.z()/Etau1,	0,	0,	0,	-ptau2.x()/Etau2,	-ptau2.y()/Etau2,	-ptau2.z()/Etau2,	0,	0,	0,	0,	0,	0,	0 };

  TMatrixD H(dimM,dimX);
  for(int i = 0; i < dimX; i++)
  {
    H(0,i) = h1[i];
    H(1,i) = h2[i];
    H(2,i) = h3[i];
    H(3,i) = h4[i];
    H(4,i) = h5[i];
    H(5,i) = h6[i];
    H(6,i) = h7[i];
    H(7,i) = h8[i];
    H(8,i) = h9[i];
    H(9,i) = h10[i];
    H(10,i) = h11[i];
    H(11,i) = h12[i];
    H(12,i) = h13[i];
    H(13,i) = h14[i];
    H(14,i) = h15[i];
    H(15,i) = h16[i];
    H(16,i) = h17[i];
    H(17,i) = h18[i];
    H(18,i) = h19[i];
    H(19,i) = h20[i];
    H(20,i) = h21[i];
    H(21,i) = h22[i];
    H(22,i) = h23[i];
    H(23,i) = h24[i];
  }
  return H;
}

TVectorD h( TVectorD x )
{
  // h(x)
  ROOT::Math::XYZPoint BV( x(0), x(1), x(2) );
  ROOT::Math::XYZVector pB( x(3), x(4), x(5) );
  double EB = x(6);
  ROOT::Math::XYZVector ptau1( x(7), x(8), x(9) );
  ROOT::Math::XYZVector pnu1( x(10), x(11), x(12) );
  ROOT::Math::XYZVector ptau2( x(13), x(14), x(15) );
  ROOT::Math::XYZVector pnu2( x(16), x(17), x(18) );
  double L1 = x(19);
  double L2 = x(20);
  double L = x(21);
  double LK = x(22);

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Enu1 = sqrt( pnu1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );

  ROOT::Math::XYZPoint PV( BV.x() - L*pB.x(), BV.y() - L*pB.y(), BV.z() - L*pB.z() );
  ROOT::Math::XYZPoint DV1( BV.x() + L1*ptau1.x(), BV.y() + L1*ptau1.y(), BV.z() + L1*ptau1.z() );
  ROOT::Math::XYZVector p3pi1( ptau1.x() - pnu1.x(), ptau1.y() - pnu1.y(), ptau1.z() - pnu1.z() );
  double E3pi1 = Etau1 - Enu1;
  ROOT::Math::XYZPoint DV2( BV.x() + L2*ptau2.x(), BV.y() + L2*ptau2.y(), BV.z() + L2*ptau2.z() );
  ROOT::Math::XYZVector p3pi2( ptau2.x() - pnu2.x(), ptau2.y() - pnu2.y(), ptau2.z() - pnu2.z() );
  double E3pi2 = Etau2 - Enu2;
  ROOT::Math::XYZPoint RP( BV.x() - LK*( pB.x() - ptau1.x() - ptau2.x() ), BV.y() - LK*( pB.y() - ptau1.y() - ptau2.y() ), BV.z() - LK*( pB.z() - ptau1.z() - ptau2.z() ) );
  ROOT::Math::XYZVector pK( pB.x() - ptau1.x() - ptau2.x(), pB.y() - ptau1.y() - ptau2.y(), pB.z() - ptau1.z() - ptau2.z() );
  double EK = EB - Etau1 - Etau2;

  TVectorD h( dimM );
  h(0) = PV.x();
  h(1) = PV.y();
  h(2) = PV.z();
  h(3) = DV1.x();
  h(4) = DV1.y();
  h(5) = DV1.z();
  h(6) = p3pi1.x();
  h(7) = p3pi1.y();
  h(8) = p3pi1.z();
  h(9) = E3pi1;
  h(10) = DV2.x();
  h(11) = DV2.y();
  h(12) = DV2.z();
  h(13) = p3pi2.x();
  h(14) = p3pi2.y();
  h(15) = p3pi2.z();
  h(16) = E3pi2;
  h(17) = RP.x();
  h(18) = RP.y();
  h(19) = RP.z();
  h(20) = pK.x();
  h(21) = pK.y();
  h(22) = pK.z();
  h(23) = EK;

  return h;
}

TVectorD x_initial_estimate( TVectorD m )
{
  // Builds an initial estimate for x based on the analytical calculations and on the known parameters in m

  ROOT::Math::XYZPoint PV( m(0), m(1), m(2) );
  ROOT::Math::XYZPoint DV1( m(3), m(4), m(5) );
  ROOT::Math::XYZVector p3pi1( m(6), m(7), m(8) );
  double E3pi1 = m(9);
  ROOT::Math::XYZPoint DV2( m(10), m(11), m(12) );
  ROOT::Math::XYZVector p3pi2( m(13), m(14), m(15) );
  double E3pi2 = m(16);
  ROOT::Math::XYZPoint RP( m(17), m(18), m(19) );
  ROOT::Math::XYZVector pK( m(20), m(21), m(22) );
  double EK = m(23);

  double m3pi1 = sqrt( pow(E3pi1,2) - p3pi1.Mag2() );
  double m3pi2 = sqrt( pow(E3pi2,2) - p3pi2.Mag2() );

  ROOT::Math::XYZPoint PV_t = makeTransformation_point( pK, RP, PV, false );
  ROOT::Math::XYZPoint DV1_t = makeTransformation_point( pK, RP, DV1, false );
  ROOT::Math::XYZVector p3pi1_t = makeTransformation_vec( pK, RP, p3pi1, false );
  ROOT::Math::XYZPoint DV2_t = makeTransformation_point( pK, RP, DV2, false );
  ROOT::Math::XYZVector p3pi2_t = makeTransformation_vec( pK, RP, p3pi2, false );
  ROOT::Math::XYZVector pK_t = makeTransformation_vec( pK, RP, pK, false );

  double a1 = DV1_t.y()/DV1_t.x();
  double a2 = DV2_t.y()/DV2_t.x();
  double b = ( PV_t.y() -a1*PV_t.x() )/( a2*PV_t.x() - PV_t.y() );
  double c = b*(DV1_t.x()/DV2_t.x());
  double d = b*(DV2_t.z() - DV1_t.z())/(DV2_t.x());
  double e = ( (1+b)*(DV1_t.z() - PV_t.z()) + d*PV_t.x() )/( (1+b)*DV1_t.x() - (1+c)*PV_t.x() );
  double f = ( PV_t.x()*sqrt(pK_t.Mag2()) )/( (1+b)*DV1_t.x() -(1+c)*PV_t.x() );
  double g = c*e + d;
  double h = f*c;
  double i = DV1_t.z() - e*DV1_t.x();
  double j =  f*DV1_t.x();

  double x1 = p3pi1_t.x() + a1*p3pi1_t.y() + e*p3pi1_t.z();
  double x2 = b*p3pi2_t.x() + a2*b*p3pi2_t.y() + g*p3pi2_t.z();

  double p1 = 1 + pow(a1,2) + pow(e,2) - pow(x1/E3pi1,2);
  double p2 = 2*e*f - ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )*(x1/pow(E3pi1,2) );
  double p3 = pow(mtau,2) + pow(f,2) - pow( ( pow(mtau,2) + pow(m3pi1,2) + 2*f*p3pi1_t.z() )/(2*E3pi1), 2);
  double q1 = pow(b,2) + pow(a2*b,2) + pow(g,2) - pow(x2/E3pi2,2);
  double q2 = 2*g*h - ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )*(x2/pow(E3pi2,2) );
  double q3 =  pow(mtau,2) + pow(h,2) - pow( ( pow(mtau,2) + pow(m3pi2,2) + 2*h*p3pi2_t.z() )/(2*E3pi2),2 );

  Double_t Ptau1x_t = (p1*q3 - p3*q1)/(p2*q1 - p1*q2);
  Double_t Ptau1y_t = a1*Ptau1x_t;
  Double_t Ptau1z_t = e*Ptau1x_t + f;
  Double_t Ptau2x_t = b*Ptau1x_t;
  Double_t Ptau2y_t = a2*b*Ptau1x_t;
  Double_t Ptau2z_t = g*Ptau1x_t + h;
  Double_t BVz_t = i - j*(1/Ptau1x_t);

  ROOT::Math::XYZVector ptau1_t( Ptau1x_t, Ptau1y_t, Ptau1z_t );
  ROOT::Math::XYZVector ptau2_t( Ptau2x_t, Ptau2y_t, Ptau2z_t );
  ROOT::Math::XYZPoint BV_t(0, 0, BVz_t);

  // transform back to LHCb frame
  ROOT::Math::XYZVector ptau1 = makeTransformation_vec( pK, RP, ptau1_t, true );
  ROOT::Math::XYZVector ptau2 = makeTransformation_vec( pK, RP, ptau2_t, true );
  ROOT::Math::XYZPoint BV = makeTransformation_point( pK, RP, BV_t, true );

  ROOT::Math::XYZVector pnu1 = ptau1 - p3pi1;
  ROOT::Math::XYZVector pnu2 = ptau2 - p3pi2;
  ROOT::Math::XYZVector pB = ptau1 + ptau2 + pK + pnu1 + pnu2;

  double Etau1 = sqrt( pow(mtau,2) + ptau1.Mag2() );
  double Etau2 = sqrt( pow(mtau,2) + ptau2.Mag2() );
  double Enu1 = sqrt( pnu1.Mag2() );
  double Enu2 = sqrt( pnu2.Mag2() );
  double EB = Etau1 + Etau2 + EK;

  double L1 = sqrt( (DV1 - BV).Mag2() )/sqrt( ptau1.Mag2() );
  double L2 = sqrt( (DV2 - BV).Mag2() )/sqrt( ptau2.Mag2() );
  double L = sqrt( (BV - PV).Mag2() )/sqrt( pB.Mag2() );
  double LK = sqrt( (RP - BV).Mag2() )/sqrt( pK.Mag2() );

  TVectorD x0(dimX);
  x0(0) = BV.x();
  x0(1) = BV.y();
  x0(2) = BV.z();
  x0(3) = pB.x();
  x0(4) = pB.y();
  x0(5) = pB.z();
  x0(6) = EB;
  x0(7) = ptau1.x();
  x0(8) = ptau1.y();
  x0(9) = ptau1.z();
  x0(10) = pnu1.x();
  x0(11) = pnu1.y();
  x0(12) = pnu1.z();
  x0(13) = ptau2.x();
  x0(14) = ptau2.y();
  x0(15) = ptau2.z();
  x0(16) = pnu2.x();
  x0(17) = pnu2.y();
  x0(18) = pnu2.z();
  x0(19) = L1;
  x0(20) = L2;
  x0(21) = L;
  x0(22) = LK;

  return x0;
}

ROOT::Math::XYZVector makeTransformation_vec(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZVector theVector, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZVector theVector_t;
  if(invFlag)
    theVector_t = myTransformation_inv * theVector;
  else 
    theVector_t = myTransformation * theVector;
  
  return theVector_t;
}

ROOT::Math::XYZPoint makeTransformation_point(ROOT::Math::XYZVector Pk, ROOT::Math::XYZPoint refPoint, ROOT::Math::XYZPoint thePoint, bool invFlag)
{
  double deltaX = refPoint.x();
  double deltaY = refPoint.y();
  double deltaZ = refPoint.z();

  ROOT::Math::Translation3D myShift(-deltaX, -deltaY, -deltaZ);

  //Using the Euler angle formalism to define the rotations
  //See https://mathworld.wolfram.com/EulerAngles.html
  //https://root.cern.ch/doc/master/classROOT_1_1Math_1_1EulerAngles.html

  //Rotation about original Z axis to bring Pk into YZ plane. If Y component is +ve, rotation is clockwise, else anti-clockwise
  double phi = -1 * TMath::ATan(Pk.x()/Pk.y());
  //Clockwise rotation about new X axis to align Z axis with Pk
  double theta = -1 * TMath::ATan(Pk.Rho() * (Pk.y()/TMath::Abs(Pk.y())) /Pk.z());

  ROOT::Math::EulerAngles myRotation(phi, theta, 0);

  //Combine Translation and EulerAngles into one single Transform 3D object

  ROOT::Math::Transform3D myTransformation     = ROOT::Math::Transform3D(myRotation) * ROOT::Math::Transform3D(myShift);
  ROOT::Math::Transform3D myTransformation_inv = myTransformation.Inverse();

  ROOT::Math::XYZPoint thePoint_t;
  if(invFlag)
    thePoint_t = myTransformation_inv * thePoint;
  else 
    thePoint_t = myTransformation * thePoint;
  
  return thePoint_t;
}
