#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <tuple>
#include <numeric>

using std::vector;
using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0,1);

double  linint( double x1, double y1, double x2, double y2, double x );
double  logint( double x1, double y1, double x2, double y2, double x );


void mass_exp (double sn_mass, int x, int y, int z, double  cell_size, int cell_no, vector<vector<vector<vector<double> > > > &iso_mass);
void print_results (int cell_no, vector<vector<vector<vector<double> > > > &iso_mass, int time, vector<vector<vector<vector<double> > > > &iso_star, vector<vector<vector<vector<double> > > > &life_star, vector<vector<vector<vector<double> > > > &age_star);
double  yields();
int i,j,k,l,m;
int Z_num; // Number of metallicities.
double  y_Z[5]; // Metallicity.
int M_num[5]; // Number of masses.
double  stellar_mass[5][24]; // List of masses for each metallicity.
double  stellar_iso[5][24][9]; // List of yields for each mass.
int elmax = 9;
double  metallicity;
const double PI  = 3.141592653589793238463;

double  powerlaw(double  power, double  lo, double  hi);
void seed();
double  salpeter(double  lo, double  hi);
int t_steps = 400;


int main() {
    t_steps++;
    int time;
    float  cell_size = 25; // Setting the size of each cell as 50 cubic parsecs.
    int cell_no = 100;
    float  time_gyr;
    int x,y,z;
    double mass_add;
    double  sfe = 0.3e-3; // 0.25; Argast // 0.02; Original
    int buffer = 5;
    double  centre = cell_no/2;
    double  mass_up = 50;
    double  star_maxmass = 40.0;
    double  star_minmass = 0.1;
    double  temp_rand;
    double gas_mass;
    int new_stars;
    double  count;
    int SN_Rate[t_steps];
    int iso;
    double total_mass;
    double distancei, distancej, distancek, rho_0, rho, radius;
    bool density_profile;

    
    density_profile = false;

    vector<vector<vector<double> > > total_stellar_mass;

    vector<vector<vector<int> > > stars;
    vector<vector<vector<vector<double> > > > star_mass;
    vector<vector<vector<vector<double> > > > life_star;
    vector<vector<vector<vector<double> > > > age_star;
    vector<vector<vector<vector<double> > > > iso_star;
    vector<vector<vector<vector<double> > > > iso_mass;

    // Setting up vector arrays...
    total_stellar_mass.resize(cell_no);
    stars.resize(cell_no);
    star_mass.resize(cell_no);
    life_star.resize(cell_no);
    age_star.resize(cell_no);
    iso_mass.resize(cell_no);
    iso_star.resize(cell_no);
    for (int i = 0; i < cell_no; ++i) {
        total_stellar_mass[i].resize(cell_no);
        stars[i].resize(cell_no);
        star_mass[i].resize(cell_no);
        life_star[i].resize(cell_no);
        age_star[i].resize(cell_no);
        iso_mass[i].resize(cell_no);
        iso_star[i].resize(cell_no);
        for (int j = 0; j < cell_no; ++j){
            total_stellar_mass[i][j].resize(cell_no);
            stars[i][j].resize(cell_no);
            star_mass[i][j].resize(cell_no);
            life_star[i][j].resize(cell_no);
            age_star[i][j].resize(cell_no);
            iso_mass[i][j].resize(cell_no);
            iso_star[i][j].resize(cell_no);
            for (int k = 0; k < cell_no; ++k) {
                star_mass[i][j][k].resize(cell_no);
                life_star[i][j][k].resize(cell_no);
                age_star[i][j][k].resize(cell_no);
                iso_mass[i][j][k].resize(cell_no);
                iso_star[i][j][k].resize(cell_no);
            }
        }
    }
    for (i=0; i<t_steps; i++){
        SN_Rate[i] = 0;
    }
    double mtot;
    double max_mass;
    int star_n;
    double SFR[t_steps];
    double  rho_c = 2.471e-2;
    rho_0 = 2e-2;
    mtot = 0;
    yields();
    for (i=0; i < cell_no; i++) {
        for (j=0; j < cell_no; j++) {
            for (k=0; k < cell_no; k++) {
                for (iso=2; iso<elmax; iso++) {
                    iso_mass[i][j][k][iso] = 1e-8;
                }
                total_stellar_mass.at(i).at(j).at(k) = 0;
                stars[i][j][k] = 0;
                for (l=0; l<cell_no; l++) {
                    star_mass.at(i).at(j).at(k).at(l) = 0;
                }
                if (density_profile == true){
                    // Finding distance from centre of volume in pc.
                    distancei = 1e-3*cell_size * abs((cell_no/2) - i);
                    distancej = 1e-3*cell_size * abs((cell_no/2) - j);
                    distancek = 1e-3*cell_size * abs((cell_no/2) - k);
                    
                    radius  = sqrt(pow(distancei,2) + pow(distancej,2) + pow(distancek,2));
                    rho = rho_0 * exp(-3*radius);
                    mtot += rho*pow(cell_size,3);
                    //rho_0 -= rho;
                    iso_mass[i][j][k][0] = 0.74 * rho * pow(cell_size,3);
                    iso_mass[i][j][k][1] = 0.26 * rho * pow(cell_size,3);
                } else {
                    iso_mass[i][j][k][0] = 0.74 * pow(10,8)/pow(cell_no,3);
                    iso_mass[i][j][k][1] = 0.26 * pow(10,8)/pow(cell_no,3);
                }
            }
        }
    }
    cout << "Total mass (if using density profile): " << mtot << endl;
    for ( time=0; time < t_steps; time++ ){
        count = 0;
        time_gyr = time/1e3;
        mass_add = 1.3447e4 * pow(time_gyr,0.4) * exp(-(time_gyr)/5.0) / (pow(cell_no, 3)); // Integrates to 10^5.

        max_mass = 0;
        total_mass = 0.0;
        // Finding the most massive cell and calculating the SFR.
        for (i=0; i < cell_no; i++) {
            for (j=0; j < cell_no; j++) {
                for (k=0; k < cell_no; k++) {
                    if ((iso_mass.at(i).at(j).at(k).at(0) + iso_mass.at(i).at(j).at(k).at(1) + iso_mass.at(i).at(j).at(k).at(8)) >= max_mass) {
                        max_mass = iso_mass.at(i).at(j).at(k).at(0) + iso_mass.at(i).at(j).at(k).at(1) + iso_mass.at(i).at(j).at(k).at(8);
                    }
                    // Infall mass_added should be primordial
//                    iso_mass[i][j][k][0] += mass_add*0.74;
//                    iso_mass[i][j][k][1] += mass_add*0.26;

                    for (l=0;l<elmax;l++){
                        total_mass += iso_mass[i][j][k][l];
                    }
                    // Either do this here, or average the density over all cells to find the average SFR and calculate the number of stars then (like Ben's)...
                    // or could also calculate the SFR for each individual cell, then use that to find the number of stars forming in each cell at each timestep (with an efficiency factor).
                    //SFR[time] += (cell_size/1000)*cell_no*2.5*pow(10,-4) * pow((mass.at(i).at(j).at(k)),1.4); // /pow(cell_size,3)
                }
            }
        }
        //SFR[time] = sfe* pow(total_mass/pow(cell_no*cell_size,3),1.4) * 1e6 * cell_size*cell_no*0.0006;//*2.5e-4 ; // 1e6 for Myr
        SFR[time] = (sfe/(cell_no*cell_no*cell_no)) * (1/pow(cell_size*cell_size*cell_size*rho_c,1.4)) * pow(total_mass,1.4) *1e6;
        if (time%10 == 0){
            cout << "SFR " << SFR[time] << endl;
            cout << "Total mass " << total_mass << endl;
        }

        new_stars = 0; // counter for use later
        // Divide by average mass of Salpeter IMF (power -2.35) to get total number of stars born per timestep.
        star_n = round(SFR[time]/0.35); // 0.35 is the average mass of from a Salpeter IMF.

        
        while (new_stars < star_n) {
            i = round(dis(gen) * (cell_no-1));
            j = round(dis(gen) * (cell_no-1));
            k = round(dis(gen) * (cell_no-1));
            // Randomly select cells to form stars until the number of stars for the timestep is reached. Scale probability of cell being chosen by mass or density of the cell.
            gas_mass = iso_mass.at(i).at(j).at(k).at(0) + iso_mass.at(i).at(j).at(k).at(1) + iso_mass.at(i).at(j).at(k).at(8); // gas mass = m(H) + m(He) + m(Z)
            if ( gas_mass >= mass_up ) {
                    // Now finding the number of stars in each cell.
                if ( dis(gen) < (gas_mass/max_mass)) {
                    stars[i][j][k]++;
                    new_stars++;
                    
                    if ( stars[i][j][k] > star_mass[i][j][k].size() ) {
                        star_mass[i][j][k].resize(stars[i][j][k]);
                        life_star[i][j][k].resize(stars[i][j][k]);
                        age_star[i][j][k].resize(stars[i][j][k]);
                   }
                    if ( iso_star.at(i).at(j).at(k).size() < (elmax+1)*(stars[i][j][k]+1) + elmax+1 ) {
                        iso_star.at(i).at(j).at(k).resize((elmax+1)*(stars[i][j][k]+1) + elmax+1);
                    }
                    if ( stars[i][j][k] > 0 ){
                        for (l=0; l<stars[i][j][k]; l++) {
                                // Check that the current array position doesn't already have a star in it.
                                // Calculating stellar lifetimes in units of Myrs.
                                if ( star_mass[i][j][k][l] == 0 ) {
                                    star_mass[i][j][k][l] = salpeter(star_minmass,star_maxmass);
                                    // Calculate mass percentage of H, He .. Z in gas, multiply percentages by mass of star formed and subtract from iso_mass...
                                    for (int m=0; m<9; m++) {
                                        if (iso_mass[i][j][k][m] > (1e-20 + star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass))) {
                                            iso_mass[i][j][k][m] -= star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                            iso_star[i][j][k][l*elmax + m] = star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                        }else{
                                            iso_mass[i][j][k][m] = 0.0;
                                            iso_star[i][j][k][l*elmax + m] = star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                        }
                                        metallicity = iso_mass[i][j][k][8];
                                    }
                                    if (metallicity==0){
                                        metallicity = 1e-10;
                                    }
                                    life_star[i][j][k][l] = round(pow(star_mass[i][j][k][l],-2.0)*11.7e3 );
                                    age_star[i][j][k][l] = time;
                                    total_stellar_mass[i][j][k] += star_mass[i][j][k][l];
                                }
                            }
                    }
                }
            }
        }
    
//        double mass_tte;
//        double mass_tta;
//        mass_tte=0;
//        mass_tta=0;
//        for (i=0; i< cell_no; i++){
//            for (j=0; j< cell_no; j++){
//                for (k=0; k<cell_no; k++){
//                    mass_tte += iso_mass[i][j][k][0] + iso_mass[i][j][k][1] + iso_mass[i][j][k][8];
//                }
//            }
//        }
        
        
        
        ///////// Checking each cell to find stars that are dying /////////
        for ( x=0; x<cell_no; x++ ){
            for ( y = 0; y<cell_no; y++ ){
                for ( z = 0; z<cell_no; z++){
                    if ( stars.at(x).at(y).at(z) > 0 ){
                        for (l = 0; l < stars.at(x).at(y).at(z); l++) {
                                if ( life_star.at(x).at(y).at(z).at(l) <= 0 && star_mass.at(x).at(y).at(z).at(l) != 0 ){ // If a star is about to die.
                                    if (star_mass.at(x).at(y).at(z).at(l) >= 10){
                                        SN_Rate[time]++;
                                    }
                                    mass_exp(star_mass.at(x).at(y).at(z).at(l),x,y,z,cell_size,cell_no, iso_mass);

                                    star_mass.at(x).at(y).at(z).at(l) = 0.0;
                                    for (m = 0; m < elmax; m++) {
                                        iso_star.at(x).at(y).at(z).at(l*elmax + m) = 0.0;
                                    }
                                    stars.at(x).at(y).at(z)--;
                                    total_stellar_mass[x][y][z] -= star_mass.at(x).at(y).at(z).at(l);

                                        
                                }

                            life_star.at(x).at(y).at(z).at(l)--;

                        }
                    }
                }
            }
        }
//        for (i=0; i< cell_no; i++){
//            for (j=0; j< cell_no; j++){
//                for (k=0; k<cell_no; k++){
//                    mass_tta += iso_mass[i][j][k][0] + iso_mass[i][j][k][1] + iso_mass[i][j][k][8];
//                }
//            }
//        }
//        cout << "Mass difference " << mass_tte - mass_tta << endl;
        
        if ( time % 1 == 0){
            cout << "Printing generation " << time << ".\n";
            print_results(cell_no, iso_mass, time, iso_star, life_star, age_star);
        }
        
    }

    ofstream results;
    results.open ("sfr_history.dat");

    for (i=0; i<t_steps; i++) {
        results << SFR[i] << endl;
    }
    results.close();
    
    results.open ("sn_rate.dat");
    
    for (i=0; i<t_steps; i++) {
        results << SN_Rate[i] << endl;
    }
    results.close();

    return 0;
}

void mass_exp (double sn_mass, int x, int y, int z, double  cell_size, int cell_no, vector<vector<vector<vector<double> > > > &iso_mass) {
    int sn_reach[343000][3] = {0};
    int sn_bound[343000][3] = {0};
    int cells_interior, boundary_cells;
    double mass_swept=0.0;
    int radius=0;
    double distancex, distancey, distancez, distance;
    float mass_interior;
    float m_rem;
    int l, m, n;
    int Z_l, Z_u; // Upper and lower Z values (closest to true metallicity).
    int M_l, M_u; // Closest masses to true stellar mass.
    double iso_exp_Zl[9];
    double iso_exp_Zu[9];
    double  iso_exp[9];
    double iso_swept[9];
    double iso_int[9];
    int bel;
    int iso,p;
    int iconst, jconst;
    cells_interior = 0;
    boundary_cells = 0;
    bool int_m, int_z;
    int_m = true;
    int_z = true;
    double mass_interior_orig;

    
    for (m=0;m<elmax;m++){
        iso_int[m] = 0;
        iso_swept[m] = 0;
    }
    // Add the iso masses from the yield tables to the cell before explosion... Then explode and redistribute the masses.
    // Finding how many cells are needed to reach the 5E4 M_sun mass pushed out from the explosion.
    
    // First, interpolate between masses, at the two nearest metallicities... Then interpolate between the two metallicities.
    for (i=0; i<Z_num; i++) {
        // Give the i and j values of the mass/metallicity locations.
        iconst = 1;
        jconst = 1;
        if (i == 0) {
            iconst = 0;
        }
        if (j == 0) {
            jconst = 0;
        }
        if (metallicity == y_Z[i]) {
            int_z = false;
            Z_l = i;
            Z_u = Z_l;
        }
        if ( y_Z[i] < metallicity ){
            Z_l = i;
        }
        if (y_Z[i] > metallicity && y_Z[i-iconst] < metallicity ) {
            Z_u = i;
        }
        for (j=0; j<M_num[i]; j++) {
            if (sn_mass == stellar_mass[i][j]) {
                int_m = false;
                M_l = i;
                M_u = M_l;
            }
            if ( stellar_mass[i][j] < sn_mass ){
                M_l = j;
            }
            if ( stellar_mass[i][j] > sn_mass && stellar_mass[i][j-jconst] < sn_mass ) {
                M_u = j;
            }
        }
    }
    if (int_m == true){
        // Find the next upper/lower mass bin, for now these can be done a metallicity of Z_l.
        // interpolate all yield values between those masses.
        for (iso=0;iso<elmax;iso++) {
            iso_exp_Zl[iso] = linint(stellar_mass[Z_l][M_l],stellar_iso[Z_l][M_l][iso],stellar_mass[Z_l][M_u],stellar_iso[Z_l][M_u][iso],sn_mass);
            if (int_z == true) {
                iso_exp_Zu[iso] = linint(stellar_mass[Z_u][M_l],stellar_iso[Z_u][M_l][iso],stellar_mass[Z_u][M_u],stellar_iso[Z_u][M_u][iso],sn_mass);
                // Interpolate between Z values.
                iso_exp[iso] = linint(y_Z[Z_l], iso_exp_Zl[iso],y_Z[Z_u], iso_exp_Zu[iso],metallicity);
            }else{
                iso_exp[iso] = iso_exp_Zl[iso];
            }
            
        }
        
    }
    
    if (int_z == false && int_m == false) {
        // Set to be the exact un-interpolated value.
        for (iso=0; iso<elmax; iso++) {
            iso_exp[iso] = stellar_iso[Z_l][M_l][iso];
        }
    }
    for (iso=0;iso<elmax; iso++) {
        if (iso_exp[iso] <=0.0) {
            iso_exp[iso] = 1e-10;
        }

    }
//    if (iso_exp[0] + iso_exp[1] + iso_exp[2] + iso_exp[3] + iso_exp[4] + iso_exp[5] + iso_exp[6] + iso_exp[7] + iso_exp[8] >= sn_mass){
//        cout << "DerP iso mass \n";
//        cout << iso_exp[0] + iso_exp[1] + iso_exp[2] + iso_exp[3] + iso_exp[4] + iso_exp[5] + iso_exp[6] + iso_exp[7] + iso_exp[8] << endl;
//        cout << "SN MASS " << sn_mass << endl;
//    }
    // Determining the mass of the remnant based upon the mass of the star (pre-SN).
    m_rem = 0;
    
    // Count outwards in a sphere from SN event until mass collected is 5E4 M_sun.
    radius = 1;

    // Dumping ejected material initially into the cell in which the star died.
    for (iso=0;iso<elmax;iso++){
        m_rem += iso_exp[m];
    }
    m_rem = sn_mass - m_rem;

    for (iso=0; iso<elmax; iso++) {
        iso_mass[x][y][z][iso] += iso_exp[iso];
    }
    double theta, phi;
    int cones;
    cones = 8;
    if (sn_mass > 10.0) { // Ensuring there aren't any AGB stars exploding as SN.
        while (mass_swept < 4.25e4) {
            mass_swept = 0.0;
            cells_interior = 0;
            boundary_cells = 0;
            for (i=x-radius; i<=x+radius; i++) {
                for (j=y-radius; j<=y+radius; j++) {
                    for (k=z-radius; k<=z+radius; k++) {
                        for (theta=0; theta<2*PI; theta+= 2*PI/cones){
                            for (phi=0; phi<PI; phi+= PI/cones) {
                                // We want the distance from the centre of a cell to be less than the radius of the sphere + cell_size/2 (approximates a sphere from cubes).
                                l=i;
                                m=j;
                                n=k;
                                // Adding if statements to allow wrapping at the edges.
                                if (i >= cell_no) {
                                    l = i - cell_no;
                                }
                                
                                if (j >= cell_no) {
                                    m = j - cell_no;
                                }
                                
                                if (k >= cell_no) {
                                    n = k - cell_no;
                                }
                                
                                if (i < 0) {
                                    l = i + cell_no;
                                }
                                
                                if (j < 0) {
                                    m = j + cell_no;
                                }
                                
                                if (k < 0) {
                                    n = k + cell_no;
                                }
                                
                                distancex = abs(i-x)*cell_size;
                                distancey = abs(j-y)*cell_size;
                                distancez = abs(k-z)*cell_size;
                                distance  = sqrt(pow(distancex,2) + pow(distancey,2) + pow(distancez,2));
                                
                                if ( (distance < radius*cell_size) && (theta >= seg*(2*PI/cones)) && (theta < (seg+1)*(2*PI/cones)) && (phi >= seg*(PI/(2*cones))) && (phi < (seg+1)*(PI/(2*cones))) ){
                                    mass_swept += iso_mass[l][m][n][0] + iso_mass[l][m][n][1] + iso_mass[l][m][n][8];
                                    // Locations of all cells interior to blast.
                                    sn_reach[cells_interior][0] = l;
                                    sn_reach[cells_interior][1] = m;
                                    sn_reach[cells_interior][2] = n;
                                    cells_interior++;
                                }
                                if (distance >= radius*cell_size  && distance < (radius+1)*cell_size && (theta >= seg*(2*PI/cones)) && (theta < (seg+1)*(2*PI/cones)) && (phi >= seg*(PI/(2*cones))) && (phi < (seg+1)*(PI/(2*cones))) ){
                                    // Locations of all cells at the boundary of the blast.
                                    sn_bound[boundary_cells][0] = l;
                                    sn_bound[boundary_cells][1] = m;
                                    sn_bound[boundary_cells][2] = n;
                                    boundary_cells++;
                                }
                                
                            }
                        }
                    }
                }
            }
            radius++;
        }

        mass_interior = 5.0; // The amount of mass to be distributed within the boundary of the explosion.
        // Now distributing that mass.

        
        mass_interior_orig = 0;
        for (i=0; i<cells_interior; i++) {
            for (m=0; m<elmax; m++) {
                iso_int[m] += iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m];
                mass_interior_orig += iso_int[m];
            }
        }

        double mass_swept;
        mass_swept = 0;
        
        for (i=0; i<cells_interior; i++) {
            // Count up all the isotope masses swept out from this region.
            for (m=0; m<elmax; m++) {
                iso_swept[m] += iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] * (1.0-(mass_interior/mass_interior_orig));
            }
            
        }
        
        for (i=0; i< cells_interior; i++) {
            for (m=0; m<elmax; m++) {
                iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] = (iso_int[m]-iso_swept[m])/cells_interior;
            }
            
        }

        for (i=0; i< boundary_cells; i++) {
            for (m=0; m<elmax; m++) {
                iso_mass[sn_bound[i][0]][sn_bound[i][1]][sn_bound[i][2]][m] += (iso_swept[m]/boundary_cells);
            }
        }

    }
    
}

void print_results( int cell_no, vector<vector<vector<vector<double> > > > &iso_mass, int time, vector<vector<vector<vector<double> > > > &iso_star, vector<vector<vector<vector<double> > > > &life_star, vector<vector<vector<vector<double> > > > &age_star) {
    int i,j,k,m;
    int slice_number = 50;
    ostringstream filename;
    ostringstream filename2;

    bool print_gas = true;
    bool print_stars = true;
    if (print_gas == true){
        filename << "gas_t" << time << ".dat";
        ofstream results;
        results.open (filename.str());

        k=slice_number;
        for (i=0; i<cell_no; i++) {
           for (j=0; j<cell_no; j++) {
                   for (m=0; m<elmax; m++) {
                       results << iso_mass[i][j][k][m] << " ";

                   }
                   results << "\n";
           }
        }
        results.close();
    }
    if (print_stars == true){
        filename2 << "stars_t" << time << ".dat";
        ofstream results;
        results.open (filename2.str());
        k=slice_number;

        for (i=0; i<cell_no; i++) {
            for (j=0; j<cell_no; j++) {
                    for (l=0; l<life_star[i][j][k].size(); l++) {
                        if (life_star[i][j][k][l] > 0.0 && iso_star[i][j][k][l*elmax+m] >= 0.0) {
                            results << time - age_star[i][j][k][l] << " "; // star age
                            for (m=0; m<elmax; m++) {
                                results << iso_star[i][j][k][l*elmax + m] << " ";
                            }
                            results << "\n";
                        }
                    }
            }
        }
        results.close();
    }
}

void seed() {
    srand(time(NULL));
}


double  powerlaw(double  power, double  lo, double  hi) {
    double  pow1,r,norm,expo,x,temp;
    pow1 = power + 1.0;
    if ( lo > hi ) {
        temp=lo;
        lo=hi;
        hi=temp;
    }
    seed();
    r = dis(gen);

    if ( power != -1.0 ) {
        norm = 1.0/(pow(hi,pow1) - pow(lo,pow1));
        expo = log10(r/norm + pow(lo,pow1))/pow1;
        x = pow(10.0,expo);
    } else {
        norm = 1.0/(log(hi) - log(lo));
        x = exp(r/norm + log(lo));
    }

    return x;
}

double  salpeter(double  lo, double  hi) {
    seed();
    return pow(((pow(hi,-1.35) - pow(lo,-1.35))*dis(gen) + pow(lo,-1.35)),(1/(-1.35)));
}

double  yields() {
    // For WW95
    /*
     Woosley_RV file format:
     
     number of metallicities
     metallicity
     number of masses
     mass H He c n o mg si Fe Z //cols 5 7 for O and Si
     */
    cout << "Reading stellar yields..." << endl;
    ifstream yield_file;
    yield_file.open("./Yields/Woosley_RV");
    
    yield_file >> Z_num;
    for (i=0; i<Z_num; i++) {
        yield_file >> y_Z[i];
        yield_file >> M_num[i];
        for (j=0; j<M_num[i]; j++) {
            yield_file >> stellar_mass[i][j];
            for (k=0; k<elmax; k++) {
                yield_file >> stellar_iso[i][j][k];
                if (stellar_iso[i][j][k] < 0.0) {
                    stellar_iso[i][j][k] = 0.0;
                }
            }
        }
    }
    yield_file.close();
    return 0;
}



/* linear interpolate/extrapolate module */
double  linint( double x1, double y1, double x2, double y2, double x )
{
    double  dydx, y_0;   /* m=slope and b=y_intercept */
    
    if( x1 == x2 ) return( ( y1 + y2 ) / 2.0 );
    if( x == x1 ) return( y1 );
    if( x == x2 ) return( y2 );
    
    return( y1 + ((x - x1)*(y2 - y1)/(x2 - x1)) );
    
}
/*******************************************************************/

/* logarithmic interpolate/extrapolate module */
double  logint( double x1, double y1, double x2, double y2, double x )
{
    double  dydx, y_0;   /* m=slope and b=y_intercept */
    
    if( x1 == x2 ) return( ( y1 + y2 ) / 2.0 );
    if( x == x1 ) return( y1 );
    if( x == x2 ) return( y2 );
    
    y1 = log10( y1 + 1.0e-8 );
    y2 = log10( y2 + 1.0e-8 );
    x1 = log10( x1 + 1.0e-8 );
    x2 = log10( x2 + 1.0e-8 );
    x = log10( x + 1.0e-8 );
    
    return( pow(10.0, (y1 + ((x - x1)*(y2 - y1)/(x2 - x1)))) - 1.0e-8 );
    
}