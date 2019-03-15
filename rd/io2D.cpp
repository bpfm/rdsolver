void open_files(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::ofstream &GENERATED_IC){
        POSITIONS.open("output/positions.txt");
        DENSITY_MAP.open("output/density.txt");
        PRESSURE_MAP.open("output/energy.txt");
        VELOCITY_MAP.open("output/momentum.txt");
        CENTRAL_COLUMN.open("output/column.txt");
        GENERATED_IC.open("output/ic.txt");
        return;
}

void close_files(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::ofstream &GENERATED_IC){
        POSITIONS.close();
        DENSITY_MAP.close();
        PRESSURE_MAP.close();
        VELOCITY_MAP.close();
        CENTRAL_COLUMN.close();
        GENERATED_IC.close();
        return;
}

void output_state(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::ofstream &GENERATED_IC, std::vector<std::vector<VERTEX> > POINTS, double T, double DT, double DX, double DY){
        int i,j;
        double TOTAL_DENSITY = 0.0;
        double TOTAL_ENERGY = 0.0;

        for(j=0;j<N_POINTS_Y;j++){
                for(i=0;i<N_POINTS_X;i++){
                        if(T == 0.0){POSITIONS << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << std::endl;}
                        GENERATED_IC << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << "\t" << POINTS[j][i].get_pressure() << "\t" << POINTS[j][i].get_x_velocity() << "\t" <<POINTS[j][i].get_y_velocity() << "\t" << POINTS[j][i].get_dx() << "\t" <<POINTS[j][i].get_dy() << std::endl;
                        DENSITY_MAP  << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << std::endl;
                        PRESSURE_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density()*POINTS[j][i].get_specific_energy() << std::endl;
                        VELOCITY_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density()*POINTS[j][i].get_x_velocity()   << "\t" << POINTS[j][i].get_mass_density()*POINTS[j][i].get_y_velocity() << std::endl;
                        if(j == N_POINTS_Y/2){CENTRAL_COLUMN << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << "\t" <<POINTS[j][i].get_pressure() << "\t" << POINTS[j][i].get_x_velocity() << "\t" << POINTS[j-1][i].get_mass_density() << "\t" << POINTS[j+1][i].get_mass_density() << std::endl;}
                        TOTAL_DENSITY += POINTS[j][i].get_mass_density()*POINTS[j][i].get_dual();
                        TOTAL_ENERGY += POINTS[j][i].get_specific_energy()*POINTS[j][i].get_dual() * POINTS[j][i].get_mass_density();
                }
        }

        std::cout << "*********************************************************" << std::endl;            // right out time and total density to terminal
        std::cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY << " \t-> total energy =\t" << TOTAL_ENERGY << "\ttime step = \t" << DT << std::endl;
        DENSITY_MAP    << " " << std::endl;
        PRESSURE_MAP   << " " << std::endl;
        VELOCITY_MAP   << " " << std::endl;
        CENTRAL_COLUMN << " " << std::endl;

        return;
}

int READ_POSITION_HEADER(std::ifstream &POSITIONS_FILE){
        std::string INFO;
        int N_POINTS;

        std::getline(POSITIONS_FILE,INFO);

        POSITIONS_FILE >> N_POINTS;

        return N_POINTS;
}

VERTEX READ_POSITION_LINE(int N_POINTS, std::ifstream &POSITIONS_FILE){
        double X,Y;
        VERTEX NEW_VERTEX;

        NEW_VERTEX.set_x(X);
        NEW_VERTEX.set_y(Y);

        return NEW_VERTEX;
}



