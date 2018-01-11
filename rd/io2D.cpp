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

void output_state(ofstream &DENSITY_MAP, ofstream &PRESSURE_MAP, ofstream &VELOCITY_MAP, ofstream &DU_FILE, vector<centre> POINTS, double T, double DT, double DX){
        int i;
        double TOTAL_DENSITY;
        vector<centre>::iterator IT_VERT;

        for(IT_VERT=POINTS.begin(),i=0;IT_VERT<POINTS.end();IT_VERT++,i++){
                        DENSITY_MAP << POINTS[i].get_x() << "\t" << POINTS[i].get_mass_density() << endl;
                        PRESSURE_MAP << POINTS[i].get_x() << "\t" << POINTS[i].get_pressure() << endl;
                        VELOCITY_MAP << POINTS[i].get_x() << "\t" << POINTS[i].get_velocity() << endl;
                        TOTAL_DENSITY += POINTS[i].get_mass_density()*DX;
        }
        cout << "*********************************************************" << endl;            // right out time and total density to terminal
        cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY  << "\ttime step = \t" << DT << endl;
        DENSITY_MAP << " " << endl;
        PRESSURE_MAP << " " << endl;
        VELOCITY_MAP << " " << endl;
        return;
}