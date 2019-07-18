void open_files(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN){
        POSITIONS.open("output/positions.txt");
        DENSITY_MAP.open("output/density.txt");
        PRESSURE_MAP.open("output/energy.txt");
        VELOCITY_MAP.open("output/momentum.txt");
        CENTRAL_COLUMN.open("output/column.txt");       
        return;
}

void close_files(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN){
        POSITIONS.close();
        DENSITY_MAP.close();
        PRESSURE_MAP.close();
        VELOCITY_MAP.close();
        CENTRAL_COLUMN.close();        
        return;
}

void output_state(std::ofstream &POSITIONS, std::ofstream &DENSITY_MAP, std::ofstream &PRESSURE_MAP, std::ofstream &VELOCITY_MAP, std::ofstream &CENTRAL_COLUMN, std::vector<VERTEX> POINTS, double T, double DT, int N_POINTS){
        int i,j;
        double TOTAL_DENSITY = 0.0;
        double TOTAL_ENERGY = 0.0;

        for(i=0;i<N_POINTS;++i){
                if(T == 0.0){POSITIONS << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << std::endl;}
                DENSITY_MAP  << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_mass_density() << std::endl;
                PRESSURE_MAP << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_mass_density()*POINTS[i].get_specific_energy() << std::endl;
                VELOCITY_MAP << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_mass_density()*POINTS[i].get_x_velocity()   << "\t" << POINTS[i].get_mass_density()*POINTS[i].get_y_velocity() << std::endl;
                if(POINTS[i].get_y() > 0.49*SIDE_LENGTH_Y and POINTS[i].get_y() < 0.51 *SIDE_LENGTH_Y){CENTRAL_COLUMN << POINTS[i].get_x() << "\t" << POINTS[i].get_y() << "\t" << POINTS[i].get_mass_density() << "\t" <<POINTS[i].get_pressure() << "\t" << POINTS[i].get_y_velocity() << std::endl;}
                TOTAL_DENSITY += POINTS[i].get_mass_density()*POINTS[i].get_dual();
                TOTAL_ENERGY += POINTS[i].get_specific_energy()*POINTS[i].get_dual() * POINTS[i].get_mass_density();
        }

        std::cout << "*********************************************************" << std::endl;            // right out time and total density to terminal
        std::cout << "time\t" << T << " \t-> total mass =\t" << TOTAL_DENSITY << " \t-> total energy =\t" << TOTAL_ENERGY << "\ttime step = \t" << DT << std::endl;
        DENSITY_MAP    << " " << std::endl;
        PRESSURE_MAP   << " " << std::endl;
        VELOCITY_MAP   << " " << std::endl;
        CENTRAL_COLUMN << " " << std::endl;

        return;
}

#ifdef QHULL_IC
int qhull_read_positions_header(std::ifstream &POSITIONS_FILE){
        std::string INFO;
        int N_POINTS;

        std::getline(POSITIONS_FILE,INFO);

        std::cout << "Details of point creation =\t" << INFO << std::endl;

        POSITIONS_FILE >> N_POINTS;

        return N_POINTS;
}

int qhull_read_triangles_header(std::ifstream &TRIANGLES_FILE){
        int N_TRIANG;

        TRIANGLES_FILE >> N_TRIANG;

        return N_TRIANG;
}

VERTEX qhull_read_positions_line(std::ifstream &POSITIONS_FILE){
        double X,Y;
        VERTEX NEW_VERTEX;

        POSITIONS_FILE >> X >> Y;

        X = SIDE_LENGTH_X*(X + 0.5);
        Y = SIDE_LENGTH_Y*(Y + 0.5);

        NEW_VERTEX = setup_vertex(X,Y);

        return NEW_VERTEX;
}

TRIANGLE qhull_read_triangles_line(std::ifstream &TRIANGLES_FILE, std::vector<VERTEX> &POINTS){
        int N_VERT,VERT0,VERT1,VERT2;
        TRIANGLE NEW_TRIANGLE;

        TRIANGLES_FILE >> N_VERT >> VERT0 >> VERT1 >> VERT2;

        // std::cout << VERT0 << "\t" << VERT1 << "\t" << VERT2 << std::endl;

        NEW_TRIANGLE.set_vertex_0(&POINTS[VERT0]);
        NEW_TRIANGLE.set_vertex_1(&POINTS[VERT1]);
        NEW_TRIANGLE.set_vertex_2(&POINTS[VERT2]);

        NEW_TRIANGLE.setup_normals();

        return NEW_TRIANGLE;
}
#endif

#ifdef CGAL_IC
int cgal_read_positions_header(std::ifstream &CGAL_FILE){
        int N_POINTS, XSHEETS, YSHEETS;
        float XLOW, YLOW, XHIGH, YHIGH;

        CGAL_FILE >> XLOW >> YLOW >> XHIGH >> YHIGH;

        // check boundaries match
        if(XLOW < 0.0 or YLOW < 0.0 or int(XHIGH) != SIDE_LENGTH_X or int(YHIGH) != SIDE_LENGTH_Y){
                std::cout << "BWARNING: Exiting on mismatching boundaries." << std::endl;
                exit(0);
        }

        CGAL_FILE >> XSHEETS >> YSHEETS;

        // check data in 1-sheet format
        if(XSHEETS !=1 or YSHEETS!=1){
                std::cout << "BWARNING: Exiting on CGAL file not in 1-sheet format." << std::endl;
                exit(0);
        }

        CGAL_FILE >> N_POINTS;

        return N_POINTS;
}

int cgal_read_triangles_header(std::ifstream &CGAL_FILE){
        int N_TRIANG;
        int EMPTY;

        CGAL_FILE >> N_TRIANG;

        return N_TRIANG;
}

VERTEX cgal_read_positions_line(std::ifstream &CGAL_FILE){
        double X,Y;
        VERTEX NEW_VERTEX;

        CGAL_FILE >> X >> Y;

        NEW_VERTEX = setup_vertex(X,Y);

        return NEW_VERTEX;
}

TRIANGLE cgal_read_triangles_line(std::ifstream &CGAL_FILE, std::vector<VERTEX> &POINTS){
        int VERT0,VERT1,VERT2;
        TRIANGLE NEW_TRIANGLE;

        CGAL_FILE >> VERT0 >> VERT1 >> VERT2;

        NEW_TRIANGLE.set_vertex_0(&POINTS[VERT0]);
        NEW_TRIANGLE.set_vertex_1(&POINTS[VERT1]);
        NEW_TRIANGLE.set_vertex_2(&POINTS[VERT2]);

        // std::cout << POINTS[VERT0].get_x() << "\t" << POINTS[VERT1].get_x() << "\t" << POINTS[VERT2].get_x() << std::endl;

        double X0,X1,X2,Y0,Y1,Y2;
        // double L1X,L1Y,L2X,L2Y,CROSS;

        // X0 = POINTS[VERT0].get_x();
        // X1 = POINTS[VERT1].get_x();
        // X2 = POINTS[VERT2].get_x();

        // Y0 = POINTS[VERT0].get_y();
        // Y1 = POINTS[VERT1].get_y();
        // Y2 = POINTS[VERT2].get_y();

        // L1X = X1 - X0;
        // L1Y = Y1 - Y0;

        // L2X = X2 - X0;
        // L2Y = Y2 - Y0;

        // CROSS = L1X*L2Y - L1Y*L2X;

        // if(CROSS < 0.0){
        //         VERT1 = VERT1 + VERT2;
        //         VERT2 = VERT1 - VERT2;
        //         VERT1 = VERT1 - VERT2;
        // }

        // check if boundary triangle

        X0 = POINTS[VERT0].get_x();
        X1 = POINTS[VERT1].get_x();
        X2 = POINTS[VERT2].get_x();

        Y0 = POINTS[VERT0].get_y();
        Y1 = POINTS[VERT1].get_y();
        Y2 = POINTS[VERT2].get_y();

        if(abs(X0 - X1) > 0.5*SIDE_LENGTH_X or abs(X0 - X2) > 0.5*SIDE_LENGTH_X or abs(X1 - X2) > 0.5*SIDE_LENGTH_X or
           abs(Y0 - Y1) > 0.5*SIDE_LENGTH_Y or abs(Y0 - Y2) > 0.5*SIDE_LENGTH_Y or abs(Y1 - Y2) > 0.5*SIDE_LENGTH_Y){
                NEW_TRIANGLE.set_boundary(1);
        }else{
                NEW_TRIANGLE.set_boundary(0);
        } 
        
        NEW_TRIANGLE.setup_normals();

        return NEW_TRIANGLE;
}

#endif



