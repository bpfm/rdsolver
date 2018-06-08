void open_files(std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN){
        DENSITY_MAP.open("output/density.txt");
        PRESSURE_MAP.open("output/pressure.txt");
        VELOCITY_MAP.open("output/velocity.txt");
        CENTRAL_COLUMN.open("output/column.txt");
        return;
}

void close_files(std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN){
        DENSITY_MAP.close();
        PRESSURE_MAP.close();
        VELOCITY_MAP.close();
        CENTRAL_COLUMN.close();
        return;
}

void output_state(std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::vector<std::vector<VERTEX> > POINTS, double T, double DT, double DX, double DY){
        int i,j;
        double TOTAL_DENSITY=0.0;

        for(j=0;j<N_POINTS_Y;j++){
                for(i=0;i<N_POINTS_X;i++){
                        DENSITY_MAP  << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << std::endl;
                        PRESSURE_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_pressure()     << std::endl;
                        VELOCITY_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_x_velocity()   << "\t" <<POINTS[j][i].get_x_velocity() << std::endl;
                        TOTAL_DENSITY += POINTS[j][i].get_mass_density()*0.5*DX*DY;
                        if(j == N_POINTS_Y/2){
                                CENTRAL_COLUMN << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_mass_density() << "\t" <<POINTS[j][i].get_pressure() << "\t" << POINTS[j][i].get_x_velocity() << std::endl;
                        }
                }
        }

        std::cout << "*********************************************************" << std::endl;            // right out time and total density to terminal
        std::cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY  << "\ttime step = \t" << DT << std::endl;
        DENSITY_MAP    << " " << std::endl;
        PRESSURE_MAP   << " " << std::endl;
        VELOCITY_MAP   << " " << std::endl;
        CENTRAL_COLUMN << " " << std::endl;

        return;
}