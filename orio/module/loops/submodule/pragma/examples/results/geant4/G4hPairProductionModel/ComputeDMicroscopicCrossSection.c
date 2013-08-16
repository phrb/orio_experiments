double G4hPairProductionModel::ComputeDMicroscopicCrossSection(
                                           double tkin,
                                           double Z,
                                           double pairEnergy)
//  differential cross section
{
  double bbbtf= 183. ;
  double bbbh = 202.4 ;
  double g1tf = 1.95e-5 ;
  double g2tf = 5.3e-5 ;
  double g1h  = 4.4e-5 ;
  double g2h  = 4.8e-5 ;

  double totalEnergy  = tkin + particleMass;
  double residEnergy  = totalEnergy - pairEnergy;
  double massratio    = particleMass/electron_mass_c2 ;
  double massratio2   = massratio*massratio ;
  double cross = 0.;

  SetElement(G4lrint(Z));

  double c3 = 0.75*sqrte*particleMass;
  if (residEnergy <= c3*z13) return cross;

  double c7 = 4.*electron_mass_c2;
  double c8 = 6.*particleMass*particleMass;
  double alf = c7/pairEnergy;
  double a3 = 1. - alf;
  if (a3 <= 0.) return cross;

  // zeta calculation
  double bbb,g1,g2;
  if( Z < 1.5 ) { bbb = bbbh ; g1 = g1h ; g2 = g2h ; }
  else          { bbb = bbbtf; g1 = g1tf; g2 = g2tf; }

  double zeta = 0;
  double zeta1 = 0.073*log(totalEnergy/(particleMass+g1*z23*totalEnergy))-0.26;
  if ( zeta1 > 0.)
  {
    double zeta2 = 0.058*log(totalEnergy/(particleMass+g2*z13*totalEnergy))-0.14;
    zeta  = zeta1/zeta2 ;
  }

  double z2 = Z*(Z+zeta);
  double screen0 = 2.*electron_mass_c2*sqrte*bbb/(z13*pairEnergy);
  double a0 = totalEnergy*residEnergy;
  double a1 = pairEnergy*pairEnergy/a0;
  double bet = 0.5*a1;
  double xi0 = 0.25*massratio2*a1;
  double del = c8/a0;

  double rta3 = sqrt(a3);
  double tmnexp = alf/(1. + rta3) + del*rta3;
  if(tmnexp >= 1.0) return cross;

  double tmn = log(tmnexp);
  double sum = 0.;

  double a4,a5,a6,a7,a8,a9,xi,xii,xi1,screen,yeu,yed,ye1,ale,cre,be,fe,ymu,ymd,ym1,alm_crm,a10,bm;
  // Gaussian integration in ln(1-ro) ( with 8 points)
  for (G4int i=0; i<8; i++)
  {
    double a4 = exp(tmn*xgi[i]);     // a4 = (1.-asymmetry)
    double a5 = a4*(2.-a4) ;
    double a6 = 1.-a5 ;
    double a7 = 1.+a6 ;
    double a9 = 3.+a6 ;
    double xi = xi0*a5 ;
    double xii = 1./xi ;
    double xi1 = 1.+xi ;
    double screen = screen0*xi1/a5 ;
    double yeu = 5.-a6+4.*bet*a7 ;
    double yed = 2.*(1.+3.*bet)*log(3.+xii)-a6-a1*(2.-a6) ;
    double ye1 = 1.+yeu/yed ;
    double ale=log(bbb/z13*sqrt(xi1*ye1)/(1.+screen*ye1)) ;
    double cre = 0.5*log(1.+2.25*z23*xi1*ye1/massratio2) ;
    double be;

    if (xi <= 1.e3) be = ((2.+a6)*(1.+bet)+xi*a9)*log(1.+xii)+(a5-bet)/xi1-a9;
    else            be = (3.-a6+a1*a7)/(2.*xi);

    double fe = (ale-cre)*be;
    if ( fe < 0.) fe = 0. ;

    double ymu = 4.+a6 +3.*bet*a7 ;
    double ymd = a7*(1.5+a1)*log(3.+xi)+1.-1.5*a6 ;
    double ym1 = 1.+ymu/ymd ;
    double alm_crm = log(bbb*massratio/(1.5*z23*(1.+screen*ym1)));
    double a10,bm;
    if ( xi >= 1.e-3)
    {
      a10 = (1.+a1)*a5 ;
      bm  = (a7*(1.+1.5*bet)-a10*xii)*log(xi1)+xi*(a5-bet)/xi1+a10;
    } else {
      bm = (5.-a6+bet*a9)*(xi/2.);
    }

    double fm = alm_crm*bm;
    if ( fm < 0.) fm = 0. ;

    sum += wgi[i]*a4*(fe+fm/massratio2);
  }

  cross = -tmn*sum*factorForCross*z2*residEnergy/(totalEnergy*pairEnergy);

  return cross;
}
