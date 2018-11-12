void open_files(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::ofstream &GENERATED_IC){
        POSITIONS.open("output/positions.txt");
        DENSITY_MAP.open("output/density.txt");
        PRESSURE_MAP.open("output/pressure.txt");
        VELOCITY_MAP.open("output/velocity.txt");
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
        double TOTAL_DENSITY=0.0;

        for(j=0;j<N_POINTS_Y;j++){
                for(i=0;i<N_POINTS_X;i++){
                        if(T == 0.0){POSITIONS << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << std::endl;}
                        GENERATED_IC << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << "\t" << POINTS[j][i].get_pressure() << "\t" << POINTS[j][i].get_x_velocity() << "\t" <<POINTS[j][i].get_y_velocity() << "\t" << POINTS[j][i].get_dx() << "\t" <<POINTS[j][i].get_dy() << std::endl;
                        DENSITY_MAP  << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_mass_density() << std::endl;
                        VELOCITY_MAP << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_y() << "\t" << POINTS[j][i].get_x_velocity()   << "\t" <<POINTS[j][i].get_y_velocity() << std::endl;
                        if(j == N_POINTS_Y/2){CENTRAL_COLUMN << POINTS[j][i].get_x() << "\t" << POINTS[j][i].get_mass_density() << "\t" <<POINTS[j][i].get_pressure() << "\t" << POINTS[j][i].get_x_velocity() << "\t" << POINTS[j-1][i].get_mass_density() << "\t" << POINTS[j+1][i].get_mass_density() << std::endl;}
                        TOTAL_DENSITY += POINTS[j][i].get_mass_density()*DX*DY;
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

VERTEX READ_IC_LINE(std::ifstream &IC_FILE){
        double X,Y,MASS_DENSITY,PRESSURE,X_VELOCITY,Y_VELOCTIY,DX,DY;
        VERTEX NEW_VERTEX;

        IC_FILE >> X >> Y >> MASS_DENSITY >> PRESSURE >> X_VELOCITY >> Y_VELOCTIY >> DX >> DY;

        //std::cout << X << "\t" << Y << "\t" << MASS_DENSITY << "\t" << PRESSURE << "\t" << X_VELOCITY << "\t" << Y_VELOCTIY << "\t" << DX << "\t" << DY << std::endl;

        NEW_VERTEX.set_x(X);
        NEW_VERTEX.set_y(Y);
        NEW_VERTEX.set_mass_density(MASS_DENSITY);                   // units kg/m^3
        NEW_VERTEX.set_x_velocity(X_VELOCITY);                       // units m/s
        NEW_VERTEX.set_y_velocity(Y_VELOCTIY);                       // units m/s
        NEW_VERTEX.set_pressure(PRESSURE);                           // units N/m^2
        NEW_VERTEX.set_dx(DX);
        NEW_VERTEX.set_dy(DY);

        NEW_VERTEX.calculate_dual();

        NEW_VERTEX.setup_specific_energy();
        NEW_VERTEX.prim_to_con();
        NEW_VERTEX.reset_du_half();
        NEW_VERTEX.reset_du();

        return NEW_VERTEX;
}