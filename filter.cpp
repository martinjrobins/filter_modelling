#include "filter.h"
#include <fstream>
void read_data_files(ComsolType &particles) {

    typedef typename ComsolType::position position;

    std::cout << "reading files..." << std::endl;
    std::ifstream pressure_file("../five layer_square/[p]_five_layer_square.txt" );
    std::ifstream vel_horz_file("../five layer_square/[vel_horz]_five_layer_square.txt" );
    std::ifstream vel_vert_file("../five layer_square/[vel_vert]_five_layer_square.txt" );
    std::string line;
    for (int i=0; i<8; ++i) {
        std::getline(pressure_file, line);
        std::getline(vel_horz_file, line);
        std::getline(vel_vert_file, line);
    }

    int i=0;
    while ( pressure_file.good() ) {
        vdouble2 pos,velocity;
        double pressure,dummy;
        std::getline(pressure_file, line);
        std::istringstream buffer(line);
        buffer >> pos[0];
        buffer >> pos[1];
        buffer >> pressure;
        buffer.clear();

        std::getline(vel_horz_file, line);
        buffer.str(line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[0];
        buffer.clear();

        std::getline(vel_vert_file, line);
        buffer.str(line);
        buffer >> dummy;
        buffer >> dummy;
        buffer >> velocity[1];
        buffer.clear();

        typename ComsolType::value_type p;
        get<position>(p) = pos;
        get<dvelocity_u>(p) = velocity[0];
        get<dvelocity_v>(p) = velocity[1];
        get<dpressure>(p) = pressure;
        if (i++ % 10 == 0) {
            particles.push_back(p);
            //std::cout << "position = "<<pos<<std::endl;
        }
    }
    std::cout << "done reading files. have "<<particles.size()<<" data points"<< std::endl;
}
