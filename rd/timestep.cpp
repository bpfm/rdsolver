void drift_update_half(int TBIN_CURRENT, int N_TRIANG, double T, double DT, std::vector<TRIANGLE> &RAND_MESH){
        int TBIN;
        for(int j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                TBIN = RAND_MESH[j].get_tbin();
                if(TBIN_CURRENT == 0 or (TBIN_CURRENT == 1 and  TBIN == 1)\
                                     or (TBIN_CURRENT == 2 and (TBIN == 2 or TBIN == 1))\
                                     or (TBIN_CURRENT == 3 and  TBIN == 1)\
                                     or (TBIN_CURRENT == 4 and (TBIN == 4 or TBIN == 2 or TBIN == 1))\
                                     or (TBIN_CURRENT == 5 and  TBIN == 1)\
                                     or (TBIN_CURRENT == 6 and (TBIN == 2 or TBIN == 1))\
                                     or (TBIN_CURRENT == 7 and  TBIN == 1)\
                                     ){
                        RAND_MESH[j].calculate_first_half(T);
                }
                // RAND_MESH[j].calculate_first_half(T);                                                 // calculate flux through TRIANGLE
                RAND_MESH[j].pass_update_half(DT);
        }
}

void jump_update_half(int TBIN_CURRENT, int N_TRIANG, double T, double DT, std::vector<TRIANGLE> &RAND_MESH){
        int TBIN;
        for(int j=0;j<N_TRIANG;++j){                                                                         // loop over all triangles in MESH
                        TBIN = RAND_MESH[j].get_tbin();
                        if(TBIN_CURRENT == 0 or (TBIN_CURRENT == 1 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 2 and (TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 3 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 4 and (TBIN == 4 or TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 5 and  TBIN == 1)\
                                             or (TBIN_CURRENT == 6 and (TBIN == 2 or TBIN == 1))\
                                             or (TBIN_CURRENT == 7 and  TBIN == 1)\
                                             ){
                                RAND_MESH[j].calculate_first_half(T);
                        }
                        if(N_TBINS == 1){
                                RAND_MESH[j].pass_update_half(DT);
                        }else if(N_TBINS == 2){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and (TBIN == 2 or TBIN == 4 or TBIN == 8)){RAND_MESH[j].pass_update_half(2.0*DT);}
                        }else if(N_TBINS == 4){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 3 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 3 and (TBIN == 4 or TBIN == 8)){RAND_MESH[j].pass_update_half(4.0*DT);}
                        }else if(N_TBINS == 8){
                                if(TBIN_CURRENT == 0 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 1 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 1 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 2 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 3 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 3 and TBIN == 4){RAND_MESH[j].pass_update_half(4.0*DT);}

                                if(TBIN_CURRENT == 4 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 5 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 5 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}

                                if(TBIN_CURRENT == 6 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}

                                if(TBIN_CURRENT == 7 and TBIN == 1){RAND_MESH[j].pass_update_half(DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 2){RAND_MESH[j].pass_update_half(2.0*DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 4){RAND_MESH[j].pass_update_half(4.0*DT);}
                                if(TBIN_CURRENT == 7 and TBIN == 8){RAND_MESH[j].pass_update_half(8.0*DT);}
                        }
                }
}