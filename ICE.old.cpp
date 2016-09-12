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
using std::vector;
using namespace std;

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0,1);

double  linint( double x1, double y1, double x2, double y2, double x );
double  logint( double x1, double y1, double x2, double y2, double x );


void mass_exp (double sn_mass, int x, int y, int z, float cell_size, int cell_no, vector<vector<vector<double> > > &mass, vector<vector<vector<vector<double> > > > &iso_mass);
void print_results (int cell_no, vector<vector<vector<double> > > &mass, vector<vector<vector<vector<double> > > > &iso_mass, int generation, vector<vector<vector<vector<double> > > > &iso_star, vector<vector<vector<vector<double> > > > &life_star);
float yields();
int i,j,k,l,m;
int Z_num; // Number of metallicities.
float y_Z[5]; // Metallicity.
int M_num[5]; // Number of masses.
float stellar_mass[5][24]; // List of masses for each metallicity.
float stellar_iso[5][24][9]; // List of yields for each mass.
int elmax = 9;
float metallicity;
const double PI  = 3.141592653589793238463;

float powerlaw(float power, float lo, float hi);
void seed();
float salpeter(float lo, float hi);
int t_steps = 10010;
int main() {
    float cell_size = 25; // Setting the size of each cell as 50 cubic parsecs.
    int cell_no = 100;
    double density = 2.83e-100;
    float time_gyr;
    int generation;
    int x,y,z;
    double mass_add;
    float nu = 0.001;
    float sfe = 0.5; // 0.25; Argast // 0.02; Original
    int buffer = 5;
    float centre = cell_no/2;
    float mass_up = 50;
    float star_maxmass = 40.0;
    float star_minmass = 0.1;
    float temp_rand;
    float gas_mass;
    float new_stars;
    float count;
    int iso;
    float total_mass;

    vector<vector<vector<double> > > mass;
    vector<vector<vector<double> > > total_stellar_mass;

    vector<vector<vector<int> > > stars;
    vector<vector<vector<vector<double> > > > star_mass;
    vector<vector<vector<vector<double> > > > life_star;
    vector<vector<vector<vector<double> > > > iso_star;
    vector<vector<vector<vector<double> > > > iso_mass;


    // Setting up vector arrays...
    mass.resize(cell_no);
    total_stellar_mass.resize(cell_no);
    stars.resize(cell_no);
    star_mass.resize(cell_no);
    life_star.resize(cell_no);
    iso_mass.resize(cell_no);
    iso_star.resize(cell_no);
    for (int i = 0; i < cell_no; ++i) {
        mass[i].resize(cell_no);
        total_stellar_mass[i].resize(cell_no);
        stars[i].resize(cell_no);
        star_mass[i].resize(cell_no);
        life_star[i].resize(cell_no);
        iso_mass[i].resize(cell_no);
        iso_star[i].resize(cell_no);
        for (int j = 0; j < cell_no; ++j){
            mass[i][j].resize(cell_no);
            total_stellar_mass[i][j].resize(cell_no);
            stars[i][j].resize(cell_no);
            star_mass[i][j].resize(cell_no);
            life_star[i][j].resize(cell_no);
            iso_mass[i][j].resize(cell_no);
            iso_star[i][j].resize(cell_no);
            for (int k = 0; k < cell_no; ++k) {
                star_mass[i][j][k].resize(cell_no);
                life_star[i][j][k].resize(cell_no);
                iso_mass[i][j][k].resize(cell_no);
                iso_star[i][j][k].resize(cell_no);
            }
        }
    }
    generation  = 0;

    float mass_count;
    float max_mass;
    int star_n;
    float SFR[t_steps];

    yields();
    for (i=0; i < cell_no; i++) {
        for (j=0; j < cell_no; j++) {
            for (k=0; k < cell_no; k++) {
                for (iso=0; iso<elmax; iso++) {
                    iso_mass[i][j][k][iso] = 1e-8;
                }
            }
        }
    }
    for ( int time=0; time < t_steps; time++ ){
        count = 0;
        time_gyr = time/1e3;
        mass_add = 1.3447e4 * pow(time_gyr,0.4) * exp(-(time_gyr)/5.0) / (pow(cell_no, 3)); // Integrates to 10^5.

        mass_count = 0;
        max_mass = 0;
        total_mass = 0;
        // Finding the most massive cell and calculating the SFR.
        for (i=0; i < cell_no; i++) {
            for (j=0; j < cell_no; j++) {
                for (k=0; k < cell_no; k++) {
                    if (mass.at(i).at(j).at(k) >= max_mass) {
                        max_mass += mass.at(i).at(j).at(k);
                    }
                    total_mass += mass.at(i).at(j).at(k);
                    // Either do this here, or average the density over all cells to find the average SFR and calculate the number of stars then (like Ben's)...
                    // or could also calculate the SFR for each individual cell, then use that to find the number of stars forming in each cell at each timestep (with an efficiency factor).
                    //SFR[time] += (cell_size/1000)*cell_no*2.5*pow(10,-4) * pow((mass.at(i).at(j).at(k)),1.4); // /pow(cell_size,3)
                }
            }
        }
        SFR[time] = (cell_size/1000)*cell_no*2.5*pow(10,-4) * pow(total_mass,1.4);
        new_stars = 0; // counter for use later
        // Divide by average mass of Salpeter IMF (power -2.35) to get total number of stars born per timestep.
        star_n = round(SFR[time]/(0.35*cell_no*cell_no*cell_no)); // 0.35 is the average mass of from a Salpeter IMF.
        //cout << "Stars " << star_n << endl;
        for (i=0; i < cell_no; i++) {
            for (j=0; j < cell_no; j++) {
                for (k=0; k < cell_no; k++) {
                    if (generation == 0){
                        mass.at(i).at(j).at(k) = pow(10,8)/pow(cell_no,3);//pow(cell_size,3) * density; // 10^8 variant copies Argast model...
                        total_stellar_mass.at(i).at(j).at(k) = 0;
                        stars[i][j][k] = 0;
                        for (l=0; l<cell_no; l++) {
                            star_mass.at(i).at(j).at(k).at(l) = 0;
                        }

                        iso_mass[i][j][k][0] = 0.74 * pow(10,8)/pow(cell_no,3);
                        iso_mass[i][j][k][1] = 0.26 * pow(10,8)/pow(cell_no,3);
                    }
                    // Infall mass_added should be primordial
//                    iso_mass[i][j][k][0] += mass_add*0.74;
//                    iso_mass[i][j][k][1] += mass_add*0.26;
                    // Randomly select cells to form stars until the number of stars for the timestep is reached. Scale probability of cell being chosen by mass or density of the cell.
//                    mass.at(i).at(j).at(k) += mass_add;
                    gas_mass = mass.at(i).at(j).at(k) -  total_stellar_mass.at(i).at(j).at(k);
                    if ( gas_mass >= mass_up ) {
                            // Now finding the number of stars in each cell.
//                            stars[i][j][k] += round(nu*pow(gas_mass,1.5)); // split_vec approach
                        if ( new_stars < star_n && dis(gen) < sfe && mass[i][j][k] >= 0.4* max_mass && (dis(gen) % star_n ) { // dis(gen) < (star_n/(cell_no*cell_no*cell_no))
                            stars[i][j][k]++;
                            new_stars ++;
                            
                        }
//                            cout << "Calculating stellar masses...\n";
                        if ( stars[i][j][k] > star_mass[i][j][k].size() ) {
                            star_mass[i][j][k].resize(stars[i][j][k]);
                            life_star[i][j][k].resize(stars[i][j][k]);
                            // Rewrite as 1D vector (accessed as 2D)... i.e. [l*elmax + m] then can just resize by adding elmax+1 spaces when l>stars[i][j][k].
                            iso_star.resize((elmax+1)*stars[i][j][k]);

                        }
                        if ( stars[i][j][k] > 0 ){
                            gas_mass = iso_mass.at(i).at(j).at(k).at(0) + iso_mass.at(i).at(j).at(k).at(1) + iso_mass.at(i).at(j).at(k).at(8); // gas mass = m(H) + m(He) + m(Z)
                        for (l=0; l<stars[i][j][k]; l++) {
//                            if ( dis(gen) < sfe) {
//                                    cout << "Stars in cell: " << stars[i][j][k] << endl;

                                try {
                                    // Check that the current array position doesn't already have a star in it.
                                    // Calculating stellar lifetimes in units of Myrs.
                                    if ( star_mass[i][j][k][l] == 0 ) {
                                        star_mass[i][j][k][l] = salpeter(star_minmass,star_maxmass);
                                        mass_count += star_mass[i][j][k][l];
                                        // Calculate mass percentage of H, He .. Z in gas, multiply percentages by mass of star formed and subtract from iso_mass...
                                        for (int m=0; m<9; m++) {
                                            if (iso_mass[i][j][k][m] > (1e-10 + star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass))) {
                                                iso_mass[i][j][k][m] -= star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                                iso_star[i][j][k][l*elmax + m] = star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                                // Correcting for subtracting Z elements twice...
                                                iso_mass[i][j][k][8] += iso_mass[i][j][k][6] + iso_mass[i][j][k][5] + iso_mass[i][j][k][4] + iso_mass[i][j][k][3] + iso_mass[i][j][k][2];
                                            }else{
                                                iso_mass[i][j][k][m] = 0.0;
                                                iso_star[i][j][k][l*elmax + m] = star_mass[i][j][k][l] * (iso_mass[i][j][k][m] / gas_mass);
                                            }
                                            metallicity = iso_mass[i][j][k][8];
                                        }
                                        if (metallicity==0){
                                            metallicity = 1e-21;
                                        }
                                        life_star[i][j][k][l] = round(pow(10,((3.79+0.24*log10(metallicity)) - (3.1+0.35*log10(metallicity))*log10(star_mass[i][j][k][l])+(0.74+0.11*log10(metallicity))*pow((log10(star_mass[i][j][k][l])),2))));
                                        total_stellar_mass[i][j][k] += star_mass[i][j][k][l];
                                    }
                                }
                                catch (exception& e) {
                                    cout << e.what() << endl;
                                    cout << "Failure" << endl;
                                }
          

//                            }else {
//                                stars[i][j][k]--;
//                            }
                            if (star_mass[i][j][k][l] >=0.01) {
                                count++;
                            }
                        }
                        }
                    }
                }
            }
        }

        ///////// Checking each cell to find stars that are dying /////////
        for ( x=0; x<cell_no; x++ ){
            for ( y = 0; y<cell_no; y++ ){
                for ( z = 0; z<cell_no; z++){
                    if ( stars.at(x).at(y).at(z) > 0 ){
                        for (l = 0; l < stars.at(x).at(y).at(z); l++) {
                            try {
                                if ( life_star.at(x).at(y).at(z).at(l) <= 0 && star_mass.at(x).at(y).at(z).at(l) != 0 ){ // If a star is about to die.
    //                                cout << "Star dying..." << star_mass[x][y][z][l] << "\n";
                                    mass_exp(star_mass.at(x).at(y).at(z).at(l),x,y,z,cell_size,cell_no, mass, iso_mass);
    //                                cout << "SN Exploded at gen " << generation << endl;
                                    star_mass.at(x).at(y).at(z).at(l) = 0.0;
                                    for (m = 0; m < elmax; m++) {
                                        iso_star.at(x).at(y).at(z).at(l*elmax + m) = 0.0;
                                    }
                                    stars.at(x).at(y).at(z)--;
                                    total_stellar_mass[x][y][z] -= star_mass.at(x).at(y).at(z).at(l);
                                }
                            }
                            catch (exception& e) {
                                cout << "Too many stars..." << endl;
                            }
                            life_star.at(x).at(y).at(z).at(l)--;

                        }
                    }
                }
            }
        }

        
        if ( generation % 10 == 0){
            cout << "Printing generation " << generation << ".\n";
            // TO DO
            // Create new array to hold values for the istope masses within the stars still alive, then pass this to the print_results routine to allow printing of stellar abundances
            print_results(cell_no, mass, iso_mass, generation, iso_star, life_star);
            cout << count << " stars formed.\n";
        }
        
        
//        cout << "SFR at gen " << generation << " is: " << (mass_count/1e6) / (2.5*2.5*2.5) << endl;
        generation++;
    }

    ofstream results;
    results.open ("sfr_history.dat");

    for (i=0; i<t_steps; i++) {
        results << SFR[i] << endl;
    }

//    total=0;
//    for (i=0; i<cell_no; i++) {
//        for (j=0; j<cell_no; j++) {
//            for (k=0; k<cell_no; k++) {
//                total += mass[i][j][k];
//            }
//        }
//    }
//    cout << "Total mass (end): " << total << endl;
//    cout << "Interior mass: " << boundmass << endl;
    return 0;
}

void mass_exp (double sn_mass, int x, int y, int z, float cell_size, int cell_no, vector<vector<vector<double> > > &mass, vector<vector<vector<vector<double> > > > &iso_mass) {

    int sn_reach[343000][3];
    int sn_bound[343000][3];
    int cells_interior, boundary_cells;
    double mass_swept=0.0;
    int radius=0;
    double distancex, distancey, distancez, distance;
    double mass_interior;
    double m_rem;
    int l, m, n;
    int Z_l, Z_u; // Upper and lower Z values (closest to true metallicity).
    int M_l, M_u; // Closest masses to true stellar mass.
    double iso_exp_Zl[9];
    double iso_exp_Zu[9];
    double iso_exp[9];
    double iso_swept[9] = {0};
    int iso;
    int iconst, jconst;
    cells_interior = 0;
    boundary_cells = 0;
    mass_interior = 0;
    bool int_m, int_z;
    int_m = true;
    int_z = true;
    mass[x][y][z] += sn_mass;
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
//                cout << "Derp? " << stellar_mass[Z_u][M_l] << " " << stellar_iso[Z_u][M_l][iso] << " " << stellar_mass[Z_u][M_u] << " " << stellar_iso[Z_u][M_u][iso] << " " << sn_mass << endl;
//                cout << "-----------------------------------------------------------------------------" << endl;
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
    // Interpolation of last isotope (i.e. iso=9 or iso=elmax) is screwing up and going to zero.
//    for (iso=0; iso<elmax; iso++) {
////        if (iso_exp[iso] <=0.0) {
////        cout << iso << " " << iso_exp[iso] << endl;
////        cout << "Yaas" << endl;
////    }
//    
//    }

//    cout << "x " << x << " y " << y << " z " << z << endl;
    // Determining the mass of the remnant based upon the mass of the star (pre-SN).
    m_rem = 3.0;
    //    double total;
    //    total=0;
    //    for (i=0; i<cell_no; i++) {
    //        for (j=0; j<cell_no; j++) {
    //            for (k=0; k<cell_no; k++) {
    //                total += mass[i][j][k];
    //            }
    //        }
    //    }
    //    cout << "Total mass (start): " << total << endl;
    
    // Count outwards in a sphere from SN event until mass collected is 5E4 M_sun.
//    cout << "Mass calculation started.\n";
    radius = 0;
    int counter;
    counter = 0;
    if (sn_mass > 10.0) { // Ensuring there aren't any AGB stars exploding as SN.
        while (mass_swept < 5e4) {
            mass_swept = 0.0;
            cells_interior = 0;
            boundary_cells = 0;

            for (i=x-radius-1; i<=x+radius+1; i++) {
                for (j=y-radius-1; j<=y+radius+1; j++) {
                    for (k=z-radius-1; k<=z+radius+1; k++) {
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
                        
                        if (distance < radius*cell_size+cell_size/2){
                            mass_swept += mass[l][m][n];
                            // Locations of all cells interior to blast.
                            sn_reach[cells_interior][0] = l;
                            sn_reach[cells_interior][1] = m;
                            sn_reach[cells_interior][2] = n;
                            cells_interior++;
                        }
                        if (distance >= radius*cell_size+cell_size/2 && distance < (radius+1)*cell_size+cell_size/2){
                            // Locations of all cells at the boundary of the blast.
                            sn_bound[boundary_cells][0] = l;
                            sn_bound[boundary_cells][1] = m;
                            sn_bound[boundary_cells][2] = n;
                            boundary_cells++;
                        }
                        if (sn_bound[boundary_cells][0]==sn_reach[cells_interior][0] && sn_bound[boundary_cells][1]==sn_reach[cells_interior][1] && sn_bound[boundary_cells][2]==sn_reach[cells_interior][2]){
                            // Ensuring boundary cells are also not used as interior cells.
                            //boundary_cells--;
                        }
                    }
                }
            }
            radius ++;
        }
    
        //    double boundmass;
    //    cout << "Mass swept out: " << mass_swept << endl;
    //    cout << "Over " << cells_interior << " cells.\n";
    //    cout << "With " << boundary_cells << " cells at the boundary \n";
        mass_interior = 5.0; // The amount of mass to be distributed within the boundary of the explosion.
        // Now distributing that mass.
        float mass_interior_orig;
        mass_interior_orig = 0;
        for (i=0; i<cells_interior; i++) {
            mass_interior_orig += mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]];
        }
        
        for (i=0; i<cells_interior; i++) {
            mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]] = mass_interior/cells_interior;
            // Count up all the isotope masses swept out from this region.
            for (m=0; m<elmax; m++) {
                iso_swept[m] += iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] * (mass_interior/mass_interior_orig);
                //if (iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m]>iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] * (mass_interior/mass_interior_orig)) {
                    iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] = (iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] * (mass_interior/mass_interior_orig));
                //cout << "BLEH " <<  (iso_mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]][m] * (mass_interior/mass_interior_orig)) << endl;
                //}
            }
            //    cout << "Mass interior " << mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]] << endl;
            //    boundmass += mass[sn_reach[i][0]][sn_reach[i][1]][sn_reach[i][2]];
            
        }
//        for (m=0; m<elmax; m++) {
//            cout << iso_swept[m] << endl;
//        }
        for (i=0; i< boundary_cells; i++) {
            mass[sn_bound[i][0]][sn_bound[i][1]][sn_bound[i][2]] += ((mass_swept-mass_interior-m_rem)/boundary_cells);
            for (m=0; m<elmax; m++) {
                iso_mass[sn_bound[i][0]][sn_bound[i][1]][sn_bound[i][2]][m] += (iso_exp[m] + iso_swept[m])/boundary_cells;
            }
            //    cout << "Mass boundary " << mass[sn_bound[i][0]][sn_bound[i][1]][sn_bound[i][2]] << endl;
        }
    }else {
        for (iso=0; iso<elmax; iso++) {
            iso_mass[x][y][z][iso] += iso_exp[iso];
        }
    }
    mass[x][y][z] += m_rem;
    for (i=0; i<343000; i++) {
        for (j=0; j<3; j++) {
            sn_reach[i][j] = 0;
            sn_bound[i][j] = 0;
        }
    }
//    return mass, iso_mass;
//    return iso_mass;
    for (i=0; i<cell_no; i++) {
        for (j=0; j<cell_no; j++) {
            for (k=0; k<cell_no; k++) {
                for (int blargle=0; blargle<elmax; blargle++) {
//                    if (iso_mass.at(i).at(j).at(k).at(blargle) < 0 ) {
//                       cout << "Mass of " << iso_mass[i][j][k][blargle] << "for iso " << blargle << endl;
//                        
//                    }
// Mn53 worth tracking fron SneIa.
//                    if (iso_mass[i][j][k][blargle] >= mass[i][j][k]) {
//                        cout << iso_exp[blargle] << endl;
//                        cout << "Isotope " << blargle << " with mass " << iso_mass[i][j][k][blargle] << endl;
//                        cout << "Gas mass " << mass[i][j][k] << endl;
//                        cout << "At coords x " << i << " y " << j << " z " << k << endl;
//                       // iso_mass[i][j][k][blargle] = 1e-10;
//
//                    }
                }
            }
        }
    }

}

void print_results( int cell_no, vector<vector<vector<double> > > &mass, vector<vector<vector<vector<double> > > > &iso_mass, int generation, vector<vector<vector<vector<double> > > > &iso_star, vector<vector<vector<vector<double> > > > &life_star) {
    int i,j,k,m;
    int slice_number = 50;
    ostringstream filename;
    ostringstream filename2;

    bool print_gas = true;
    bool print_stars = true;
    if (print_gas == true){
        filename << "gas_t" << generation << ".dat";
      //  filename += generation.str();
        ofstream results;
        results.open (filename.str());
        //    results << "X   Y   Z   Mass\n"; // with some yields now...
        // TO DO
        // print out stellar abundances, and change code so that just a single slice is written out (for slice plotting), but I suppose all stars should be written out for completeness... Thoughts...
        k=slice_number;
        for (i=0; i<cell_no; i++) {
           for (j=0; j<cell_no; j++) {
    //           for (k=0; k<cell_no; k++) {
                   results << mass[i][j][k] << " ";
                   for (m=0; m<elmax; m++) {
                       results << iso_mass[i][j][k][m] << " ";
    //                   if (iso_mass[i][j][k][m] <= 0) {
    //                       cout << "Yaasss" << endl;
    //                   }
                   }
                   results << "\n";
    //           }
           }
        }
        results.close();
    }
    if (print_stars == true){
        filename2 << "stars_t" << generation << ".dat";
        //  filename += generation.str();
        ofstream results;
        results.open (filename2.str());
        //    results << "star_age << isotopes   \n"; // with some yields now...
        // TO DO
        // print out stellar abundances, and change code so that just a single slice is written out (for slice plotting), but I suppose all stars should be written out for completeness... Thoughts...
        k=slice_number;

        for (i=0; i<cell_no; i++) {
            for (j=0; j<cell_no; j++) {
//                for (k=0; k<cell_no; k++) {
                    for (l=0; l<life_star[i][j][k].size(); l++) {
                        if (life_star[i][j][k][l] >0.0 && iso_star[i][j][k][l*elmax+m] >= 0.0) {
                         //   results << (generation - life_star[i][j][k][l]) << " "; // star age
                            for (m=0; m<elmax; m++) {
                                results << iso_star[i][j][k][l*elmax + m] << " ";
                            }
                            results << "\n";
                        }
                    }
                //}
            }
        }
        results.close();
    }
}

void seed() {
    srand(time(NULL));
}


float powerlaw(float power, float lo, float hi) {
    float pow1,r,norm,expo,x,temp;
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

float salpeter(float lo, float hi) {
    seed();
    return pow(((pow(hi,-1.35) - pow(lo,-1.35))*dis(gen) + pow(lo,-1.35)),(1/(-1.35)));
}

float yields() {
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
