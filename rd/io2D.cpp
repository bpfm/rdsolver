using namespace std;

void open_files(ofstream &DENSITY_MAP, ofstream &PRESSURE_MAP, ofstream &VELOCITY_MAP, ofstream &DU_FILE){
        DENSITY_MAP.open("density.txt");
        PRESSURE_MAP.open("pressure.txt");
        VELOCITY_MAP.open("velocity.txt");
        DU_FILE.open("mass_flux.txt");
        return;
}

void close_files(ofstream &DENSITY_MAP, ofstream &PRESSURE_MAP, ofstream &VELOCITY_MAP, ofstream &DU_FILE){
        DENSITY_MAP.close();
        PRESSURE_MAP.close();
        VELOCITY_MAP.close();
        DU_FILE.close();
        return;
}

void output_state(ofstream &DENSITY_MAP, ofstream &PRESSURE_MAP, ofstream &VELOCITY_MAP, ofstream &DU_FILE, vector<vector<VERTEX> > POINTS, double T, double DT, double DX, double DY){
        int i,j;
        double TOTAL_DENSITY;

        for(j=0;j<N_POINTS;j++){
                for(i=0;i<N_POINTS;i++){
                        DENSITY_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << endl;
                        PRESSURE_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" <<POINTS[j][i].get_pressure() << endl;
                        VELOCITY_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" <<POINTS[j][i].get_x_velocity() << "\t" <<POINTS[j][i].get_x_velocity() << endl;
                        TOTAL_DENSITY += POINTS[j][i].get_mass_density()*0.5*DX*DY;
                        //cout << "density =\t" << POINTS[j][i].get_mass_density() << "\tDX =\t" << DX << endl;
                }
        }

        cout << "*********************************************************" << endl;            // right out time and total density to terminal
        cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY  << "\ttime step = \t" << DT << endl;
        DENSITY_MAP << " " << endl;
        PRESSURE_MAP << " " << endl;
        VELOCITY_MAP << " " << endl;

        return;
}