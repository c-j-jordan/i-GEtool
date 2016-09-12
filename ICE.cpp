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
#include <algorithm>
#include <thread>

//using std::vector;
using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0,1);

double  linint( double x1, double y1, double x2, double y2, double x );
double  logint( double x1, double y1, double x2, double y2, double x );


void mass_exp (double sn_mass, int x, int y, int z, double  cell_size, int cell_no, vector<vector<vector<vector<double> > > > &iso_mass, vector<double> iso_orig);
void print_results (int cell_no, vector<vector<vector<vector<double> > > > iso_mass, int time, vector<vector<double> > iso_star, vector<double> life_star, vector<double> age_star);
double  yields();
int i,j,k,l,m;
int Z_num; // Number of metallicities.
vector<double> y_Z;
vector<int> M_num; // Number of masses.
vector<vector<double> > stellar_mass; // List of masses for each metallicity.
vector<vector<double> > m_rem; // List of remnant masses for each metallicity.
vector<vector<vector<double> > > stellar_iso; // List of yields for each mass.
vector<vector<vector<double> > > temp_stellar_iso; // List of yields for each mass.

vector<vector<vector<double> > > stellar_iso_orig;
vector<vector<vector<double> > > protno;
vector<vector<vector<double> > > massno;

int elmax = 9;
int isomax = 280;

double  metallicity;
const double PI  = 3.141592653589793238463;
int delta_t;
double  powerlaw(double  power, double  lo, double  hi);
void seed();
double  salpeter(double  lo, double  hi);
int t_steps = 10000;
int yield_choice = 2;


int main() {
    delta_t=1;
    t_steps++;
    double temp_mass;
    int time;
    float  cell_size = 25; // Setting the size of each cell as 50 cubic parsecs.
    int cell_no = 100;
    float  time_gyr;
    int x,y,z;
    double mass_add;
    double  sfe = 0.45e-4; // 0.25; Argast // 0.02; Original
    double  mass_up = 50.0;
    double  star_maxmass = 50.0;
    double  star_minmass = 0.1;
    double  temp_rand;
    double gas_mass;
    int new_stars;
    double  count;
    int SN_Rate[t_steps/delta_t];
    int iso,m;
    double total_mass;
    double distancei, distancej, distancek, rho_0, rho, radius;
    bool density_profile;

    
    density_profile = false;

    vector <double> newColumn;
    vector <int> newStar;

    vector<vector<int> > stars; // Stores the cartesian coordinates of all stars alive.
    vector<double> star_mass;
    vector<double>  life_star;
    vector<double> age_star;
    vector<vector<double> > iso_star;
    vector<vector<vector<vector<double> > > > iso_mass;

    // Setting up vector arrays...
    iso_mass.resize(cell_no);
    for ( i = 0; i < cell_no; ++i) {
        iso_mass[i].resize(cell_no);
        for ( j = 0; j < cell_no; ++j){
            iso_mass[i][j].resize(cell_no);
            for ( k = 0; k < cell_no; ++k) {
                iso_mass[i][j][k].resize(elmax);
            }
        }
    }

    int timer;
    timer=0;
    for (i=0; i<t_steps/delta_t; ++i){
        SN_Rate[timer] = 0;
        timer++;
    }
    double mtot;
    double max_mass;
    int star_n;
    double SFR[t_steps/delta_t];
    double  rho_c = 2.471e-2;
    rho_0 = 2e-2;
    mtot = 0;
    yields();
    for (i=0; i < cell_no; ++i) {
        for (j=0; j < cell_no; ++j) {
            for (k=0; k < cell_no; ++k) {
                for (iso=2; iso<elmax; ++iso) {
                    iso_mass[i][j][k][iso] = 1e-8;

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
                    iso_mass[i][j][k][0] = 6.0*0.74 * pow(10,7)/pow(cell_no,3);
                    iso_mass[i][j][k][1] = 6.0*0.26 * pow(10,7)/pow(cell_no,3);

                }
            }
        }
    }
    int old_star_size;
    cout << "Total mass (if using density profile): " << mtot << endl;
    for ( time=0; time < t_steps; time+=delta_t ){
        count = 0;
        time_gyr = time/1e3;
        //mass_add = 2.3257e7 * pow(time_gyr,0.25) * exp(-(time_gyr)/2.0) / (pow(cell_no, 3)); // Integrates to 5e7. Coeffecient 'a' solved used WolframAlpha for simplicity.
        mass_add = ((1.3447e1* pow(time_gyr,0.25) * exp(-1.0*time_gyr/1.0)) / pow(cell_no, 3)) * delta_t; // Integrates to 10^5.
        //cout << "Mass added " << mass_add << endl;
        // Infall decline timescale (tau) = 2 Gyr
        // Age of system (t_end) = 13 Gy
        // time of maximal infall (t_max) = 0.5 Gyr
        // Total infalling mass = 5e7 M_sun
        
        max_mass = 0.0;
        total_mass = 0.0;
        // Finding the most massive cell and calculating the SFR.
        //cout << time << endl;
        
        for (i=0; i < cell_no; ++i) {
            for (j=0; j < cell_no; ++j) {
                for (k=0; k < cell_no; ++k) {
                    // Infall mass_added should be primordial
                    iso_mass[i][j][k][0] += mass_add*0.74;
                    iso_mass[i][j][k][1] += mass_add*0.26;
                    temp_mass = accumulate(iso_mass[i][j][k].begin(), iso_mass[i][j][k].end(),0.0);
                    if (temp_mass >= max_mass) {
                        max_mass = temp_mass;
                    }
                    total_mass += temp_mass;
                }
            }
        }

        SFR[time] = (sfe/(cell_no*cell_no*cell_no)) * (1/pow(cell_size*cell_size*cell_size*rho_c,1.4)) * pow(total_mass,1.4) *1e6*delta_t;
        if (total_mass/(cell_no*cell_no*cell_no) < mass_up){
            SFR[time] = 0.0;
        }

        if (time%100 == 0){
            cout << "SFR " << SFR[time]/delta_t << endl;
            cout << "Total mass " << total_mass << endl;
            cout << "Infalling mass (per Myr) " << mass_add * cell_no * cell_no * cell_no << endl;
            if ( time > 1 ){
                cout << "SN Rate " << SN_Rate[time-1]/delta_t << endl;
            }
        }

        new_stars = 0; // counter for use later
        // Divide by average mass of Salpeter IMF (power -2.35) to get total number of stars born per timestep.
        star_n = round(SFR[time]); // 0.35 is the average mass of from a Salpeter IMF.
        old_star_size = stars.size();
        stars.reserve(old_star_size + new_stars);
        
        while (new_stars < star_n) {
            // Check to see if the lifetime of the star is less than the timestep size...
            i = round(dis(gen) * (cell_no-1));
            j = round(dis(gen) * (cell_no-1));
            k = round(dis(gen) * (cell_no-1));
            // Randomly select cells to form stars until the number of stars for the timestep is reached. Scale probability of cell being chosen by mass or density of the cell.
            gas_mass = accumulate(iso_mass[i][j][k].begin(), iso_mass[i][j][k].end(), 0.0); // gas mass = m(H) + m(He) + m(Z)
            if ( gas_mass >= mass_up ) {
                    // Now finding the number of stars in each cell.
                if ( dis(gen) < (pow(gas_mass,1.4)/max_mass)) {
                    stars.emplace_back(newStar);
                    stars[stars.size()-1].emplace_back(i);
                    stars[stars.size()-1].emplace_back(j);
                    stars[stars.size()-1].emplace_back(k);

                    new_stars++;
                    // Calculating stellar lifetimes in units of Myrs.
                    star_mass.emplace_back(salpeter(star_minmass,star_maxmass));
                    // Calculate mass percentage of H, He .. Z in gas, multiply percentages by mass of star formed and subtract from iso_mass...
                    for ( m=0; m<elmax; ++m) {
                        if (iso_mass[i][j][k][m] > (1e-20 + star_mass[stars.size()-1] * (iso_mass[i][j][k][m] / gas_mass))) {
                            iso_mass[i][j][k][m] -= star_mass[stars.size()-1] * (iso_mass[i][j][k][m] / gas_mass);
                            iso_star.emplace_back(newColumn);
                            iso_star[stars.size()-1].emplace_back(star_mass[stars.size()-1] * (iso_mass[i][j][k][m] / gas_mass));
                        }else{
                            iso_mass[i][j][k][m] = 0.0;
                            iso_star.emplace_back(newColumn);
                            iso_star[stars.size()-1].emplace_back(star_mass[stars.size()-1] * (iso_mass[i][j][k][m] / gas_mass));
                        }

                    }
                    metallicity = accumulate(iso_star[stars.size()-1].begin(), iso_star[stars.size()-1].end(),0.0) - (iso_star[stars.size()-1][0] + iso_star[stars.size()-1][1]);
                    //life_star.emplace_back(round(pow(star_mass[stars.size()-1],-2.0)*11.7e3 ));
                    life_star.emplace_back(round(pow(10, (3.79 + 0.24 * metallicity) - (3.10 + 0.35 * metallicity) * log10(star_mass[stars.size()-1]) + (0.74 + 0.11 * metallicity) * pow(log10(star_mass[stars.size()-1]),2) ) ) );
                    age_star.emplace_back(time);
                }
            }
        }


        vector<double> iso_orig;
        iso_orig.resize(elmax);
        vector<int> to_delete;
        to_delete.clear();

        ///////// Checking each cell to find stars that are dying /////////
        for(int death=0;death<stars.size();++death){
            // Check if next star is > 10 cells away and create new thread for that star... Should be fine to modify the same vector in different places over multiple threads.
            // Change looping variable from i to anything else...
            if ( life_star[death] <= 0 && star_mass[death] != 0.0 ){ // If a star is about to die.
                if (star_mass[death] >= 0.1){
                    SN_Rate[time]++;
                }
                if (star_mass[death] <=30.0){
                    x = stars[death][0];
                    y = stars[death][1];
                    z = stars[death][2];
                    for (m = 0; m < elmax; ++m){
                        iso_orig[m] = iso_star[death][m];
                    }

                    mass_exp(star_mass[death],x,y,z,cell_size,cell_no, iso_mass, iso_orig);


                }
                to_delete.emplace_back(death);
            }
                life_star[death]-=delta_t;
            
        }

        for (i=0;i<to_delete.size(); i++ ){
            // Erase all info of dead stars from arrays.
            iso_star.erase(iso_star.begin() + to_delete[i]);
            stars.erase(stars.begin() + to_delete[i]);
            star_mass.erase(star_mass.begin() + to_delete[i]);
            life_star.erase(life_star.begin() + to_delete[i]);
            age_star.erase(age_star.begin() + to_delete[i]);
        }
        if ( time % 100 == 0){
            cout << "Printing generation " << time << ".\n";
            print_results(cell_no, iso_mass, time, iso_star, life_star, age_star);
        }
    }

    ofstream results;
    results.open ("sfr_history.dat");

    timer=0;
    for (i=0; i<t_steps; i+= delta_t) {
        results << SFR[timer]/delta_t << endl;
        timer++;
    }
    results.close();
    
    results.open ("sn_rate.dat");
    timer=0;
    for (i=0; i<t_steps; i+=delta_t) {
        results << SN_Rate[timer]/delta_t << endl;
        timer++;
    }
    results.close();

    return 0;
}

void mass_exp (double sn_mass, int x, int y, int z, double  cell_size, int cell_no, vector<vector<vector<vector<double> > > > &iso_mass, vector<double> iso_orig) {

    vector<int> sn_reach;
    vector<int> sn_bound;
    int cells_interior, boundary_cells;
    double mass_swept;
    int radius=0;
    double distancex, distancey, distancez, distance;
    double mass_interior, mass_end;
    double m_rem;
    int l, m, n;
    int Z_l, Z_u; // Upper and lower Z values (closest to true metallicity).
    int M_l, M_u; // Closest masses to true stellar mass.
    vector<double> iso_exp_Zl;
    vector<double> iso_exp_Zu;
    vector<double>  iso_exp;
    vector<double> iso_swept;
    vector<double> iso_int;
    int iso;
    int iconst, jconst;
    bool int_m, int_z;
    int_m = true;
    int_z = true;
    double orig_mass;
    iso_exp_Zl.resize(elmax);
    iso_exp_Zu.resize(elmax);
    iso_exp.resize(elmax);
    iso_swept.resize(elmax);
    iso_int.resize(elmax);
    
    orig_mass = sn_mass;
    
    for (iso=0;iso<elmax;++iso){
        iso_int[iso] = 0;
        iso_swept[iso] = 0;
    }
    if (yield_choice == 1){
        metallicity = iso_orig[8];
    }
    if (yield_choice == 2){
        metallicity = accumulate(iso_orig.begin(), iso_orig.end(),0.0) - (iso_orig[0] + iso_orig[1]);
        // Remnant mass for stars > 25 increases exponentially allowing yields for said stars to be set to 25... As per Benoit's suggestion.
        if (sn_mass > 25.0){
            sn_mass = 25.0;
        }
    }
    // First, interpolate between masses, at the two nearest metallicities... Then interpolate between the two metallicities.
    for (i=0; i<Z_num; ++i) {
        // Give the i and j values of the mass/metallicity locations.
        iconst = i-1;
        if (i == 0) {
            iconst = 0;
        }

        if (metallicity == y_Z[i]) {
            int_z = false;
            Z_l = i;
            Z_u = Z_l;
        }
        if ( y_Z[i] < metallicity ){
            Z_l = i;
        }
        if (y_Z[i] > metallicity && y_Z[iconst] < metallicity ) {
            Z_u = i;
        }
        if (metallicity < y_Z[0]){
            Z_u = 0;
            Z_l = 0;
        }
        for (j=0; j<M_num[i]; ++j) {
            jconst = j-1;
            if (j == 0) {
                jconst = 0;
            }
            if (sn_mass == stellar_mass[i][j]) {
                int_m = false;
                M_l = j;
                M_u = M_l;
            }
            if ( stellar_mass[i][j] < sn_mass ){
                M_l = j;
            }

            if ( stellar_mass[i][j] > sn_mass && stellar_mass[i][jconst] < sn_mass ) {
                M_u = j;
            }
            if ( sn_mass < stellar_mass[i][0] ){
                M_u = 0;
                M_l = 0;
            }
        }
    }
    // Allowing extrapolation if stellar mass/metallicity is beyond the range of the yield tables.
    if (Z_l == Z_u && Z_l == 0){
        Z_u = 1;
    }
    if (Z_l == Z_u && Z_l == Z_num-1){
        Z_l = Z_num - 2;
    }
    if (M_l == M_u && M_l == 0){
        M_u = 1;
    }
    if (M_l == M_u && M_l == M_num[0]-1){
        M_l = M_num[0]-2;
    }
    if (int_m == true or int_z == true){
        // Find the next upper/lower mass bin, for now these can be done a metallicity of Z_l.
        // interpolate all yield values between those masses.
        for (iso=0;iso<elmax;++iso) {

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
        for (iso=0; iso<elmax; ++iso) {
            iso_exp[iso] = stellar_iso[Z_l][M_l][iso];
        }
    }
    for (iso=0;iso<elmax; ++iso) {
        if (iso_exp[iso] <1e-20) {
            iso_exp[iso] = 1e-20;
        }
    }
//    cout << "SN Mass " << orig_mass << endl;
//    cout << "Iso ejected " << accumulate(iso_exp.begin(), iso_exp.end(), 0.0)* (orig_mass/sn_mass) << endl;
    // Determining the mass of the remnant based upon the mass of the star (pre-SN).
    
    // Count outwards in a sphere from SN event until mass collected is 5E4 M_sun.
    radius = 1;

    // Dumping ejected material initially into the cell in which the star died.
    for (iso=0; iso<elmax; ++iso) {
        iso_mass[x][y][z][iso] += iso_exp[iso] * (orig_mass/sn_mass);
    }
    mass_interior = 0.0;
    mass_swept = 0.0;
    if (sn_mass > 10.0 or dis(gen)<1e-3) { // Ensuring there aren't any AGB stars exploding as SN (also now with 1a).
        radius=0;
        while (mass_swept < (4e4)) {
            mass_swept = 0.0;
            mass_interior = 0.0;
            cells_interior = 0;
            boundary_cells = 0;
            sn_reach.clear();
            sn_bound.clear();
            for (i=x-radius; i<=x+radius; ++i) {
                for (j=y-radius; j<=y+radius; ++j) {
                    for (k=z-radius; k<=z+radius; ++k) {
                        // Loop now over segments, just need an equation showing which cells in the sphere belong to which segment.
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

                        if ( distance < radius*cell_size ){
                            // Locations of all cells interior to blast.
                            for (iso=0; iso<elmax; iso++){
                                mass_swept += iso_mass[l][m][n][iso];
                                mass_interior += iso_mass[l][m][n][iso];
                            }

                            sn_reach.emplace_back(l);
                            sn_reach.emplace_back(m);
                            sn_reach.emplace_back(n);

                            cells_interior++;
                        }
                        if (distance >= radius*cell_size  && distance < (radius+1)*cell_size ){
                            // Locations of all cells at the boundary of the blast.
                            for (iso=0; iso<elmax; iso++){
                                mass_swept += iso_mass[l][m][n][iso];
                            }
                            sn_bound.emplace_back(l);
                            sn_bound.emplace_back(m);
                            sn_bound.emplace_back(n);
                            
                            boundary_cells++;
                        }
                    }
                }
            }
            radius++; // radius of largest segment, or radius of each segment in an array
        }
        //cout << mass_swept << endl;
        mass_end = 5.0; // The amount of mass to be distributed within the boundary of the explosion.
        // Now distributing that mass.
        
        for (i=0; i<cells_interior; ++i) {
            for (iso=0; iso<elmax; ++iso) {
                iso_int[iso] += iso_mass[sn_reach[3*i]][sn_reach[3*i+1]][sn_reach[3*i+2]][iso];
            }
        }
        for (i=0; i<cells_interior; ++i) {
            // Count up all the isotope masses swept out from this region.
            for (iso=0; iso<elmax; ++iso) {
                iso_swept[iso] += (iso_mass[sn_reach[3*i]][sn_reach[3*i+1]][sn_reach[3*i+2]][iso] * (1.0-(mass_end/mass_interior)));
            }
            
        }

        for (i=0; i<cells_interior; ++i) {
            for (iso=0; iso<elmax; ++iso) {
                iso_mass[sn_reach[3*i]][sn_reach[3*i+1]][sn_reach[3*i+2]][iso] = (iso_int[iso]-iso_swept[iso])/cells_interior;
            }
        }
        for (i=0; i<boundary_cells; ++i) {
            for (iso=0; iso<elmax; ++iso) {
                iso_mass[sn_bound[3*i]][sn_bound[3*i+1]][sn_bound[3*i+2]][iso] += (iso_swept[iso]/boundary_cells);
            }
        }

    }

}

void print_results( int cell_no, vector<vector<vector<vector<double> > > > iso_mass, int time, vector<vector<double> > iso_star, vector<double> life_star, vector<double> age_star) {
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
        for (i=0; i<cell_no; ++i) {
           for (j=0; j<cell_no; ++j) {
                   for (m=0; m<elmax; ++m) {
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

        for (i=0; i<age_star.size(); ++i) {
            if (life_star[i] > 0.0 && dis(gen) < 1) { // Only printing 10 percent of the stars.
                results << time - age_star[i] << " "; // star age
                for (m=0; m<elmax; ++m) {
                    results << iso_star[i][m] << " ";
                }
                results << "\n";
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
    // Define yield choice globally
    // also define
    // and m_rem as 2d array
    // and stellar lifetime as 2d array global.
    vector<int> elements_inuse = {0,1,4,6,10,12,25};
    double metals;
    char blank[40];
    int blah;
    cout << "Reading stellar yields..." << endl;
    ifstream yield_file;
    switch ( yield_choice ){
        case 1:
            // For WW95
            /*
             Woosley_RV file format:
             
             number of metallicities = 5
             metallicity
             number of masses = 24
             mass H He c n o mg si Fe Z //cols 5 7 for O and Si
             elmax = 9
             */
            elmax = 9;
            y_Z.resize(5);
            M_num.resize(5);
            stellar_mass.resize(5);
            stellar_iso.resize(5);
            for (i = 0; i<5; ++i){
                stellar_mass[i].resize(24);
                stellar_iso[i].resize(24);
                for (j = 0; j<24; ++j){
                    stellar_iso[i][j].resize(elmax);
                }
            }
            yield_file.open("./Yields/Woosley_RV");
            yield_file >> Z_num;
            for (i=0; i<Z_num; ++i) {
                yield_file >> y_Z[i];
                yield_file >> M_num[i];
                for (j=0; j<M_num[i]; ++j) {
                    yield_file >> stellar_mass[i][j];
                    for (k=0; k<elmax; ++k) {
                        yield_file >> stellar_iso[i][j][k];
                        if (stellar_iso[i][j][k] < 0.0) {
                            stellar_iso[i][j][k] = 0.0;
                        }
                    }
                }
            }
            yield_file.close();
            break;
        case 2:
            /* 
             NUGRID Yields as of 18/07/2016
             elmax = 280
             81 elements
             Set stars of mass greater than 25 solar masses to not explode as SN... Just include stellar winds.
            */
            cout << "Using NuGrid yield tables." << endl;
            stellar_mass.resize(5);
            stellar_iso.resize(5);
            temp_stellar_iso.resize(5);
            stellar_iso_orig.resize(5);
            protno.resize(5);
            massno.resize(5);
            for (i = 0; i<5; ++i){
                stellar_mass[i].resize(12);
                stellar_iso[i].resize(12);
                temp_stellar_iso[i].resize(12);
                stellar_iso_orig[i].resize(12);
                protno[i].resize(12);
                massno[i].resize(12);
                for (j = 0; j<12; ++j){
                    temp_stellar_iso[i][j].resize(81);
                    for (k = 0; k< 81; ++k){
                        temp_stellar_iso[i][j][k] = 0.0;
                    }
                    stellar_iso_orig[i][j].resize(280);
                    protno[i][j].resize(280);
                    massno[i][j].resize(280);
                }
            }
            
            
            yield_file.open("./Yields/isotope_yield_table_MESA_only_fryer12_delay.txt");
            yield_file >> blank >> blank >> blank >> blank >> Z_num;
            yield_file >> blank >> blank >> blank >> blank >> blank;
            y_Z.resize(Z_num);
            M_num.resize(Z_num);
            //define & m_int...
            //read values
            elmax = 0;
            int kconst;
            double temp_iso;
            for( i=Z_num-1; i >=0; --i) {
                M_num[i] = 12;
                for( j=0; j < M_num[i]; ++j) {
                    yield_file >> blank >> blank >> blank >> stellar_mass[i][j] >> blank >> y_Z[i];
//                    yield_file >> blank >> blank >> stellar_lifetime[i][j]; // Use this line if reading in stellar lifetimes from yield tables.
                    yield_file >> blank >> blank >> blank;
//                    yield_file >> blank >> blank >> m_rem[i][j]; // Use this line if you care about remnant masses.
                    yield_file >> blank >> blank >> blank;

                    yield_file >> blank >> blank >> blank >> blank >> blank;
                    elmax = 0;
                    for( k=0;k<isomax; ++k){
                        kconst = 1;
                        if (k==0){
                            kconst = 0;
                        }
                        temp_iso = 0.0;

                        yield_file >> blank >> temp_iso >> stellar_iso_orig[i][j][k] >> protno[i][j][k] >> massno[i][j][k];
                        if (protno[i][j][k] != protno[i][j][k-kconst]){
                            elmax++;
                        }
                        temp_stellar_iso[i][j][elmax] += temp_iso;

                    }
                }
            }
            elmax++;
            cout << "Total elements " << elmax << endl;
            blah = 0;
            // Selecting the elements to be used in GCE.
//            elements_inuse.emplace_back(0);
            for (i=0;i<Z_num; i++){
                for (j=0; j<M_num[i]; j++){
                    metals = 0.0;
                    for (k=0; k<elmax; k++){
                        if (find(elements_inuse.begin(), elements_inuse.end(), k) == elements_inuse.end()){
                            metals += temp_stellar_iso[i][j][k];
                        }else{
                            stellar_iso[i][j].emplace_back(temp_stellar_iso[i][j][k]);
                            if( i==0 && j==0) blah++;
                        }
                    }
                    stellar_iso[i][j].emplace_back(metals);
                }
            }
            blah++;
            temp_stellar_iso.clear();
            elmax = blah;
            cout << "Elements in use " << blah << endl;
            yield_file.close();
            break;
        default:
            cout << "Error in reading yields.\n" ;
            break;
    }
    return 0.0;
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