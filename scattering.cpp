//#pragma GCC optimize("0fast") //TODO GCC
#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <stdio.h>
#include <time.h>

using namespace std;

float const impulse_coefficient = 0.000126;
int const cell_size = 500;


struct vector_3D
{
    float x;
    float y;
    float z;
    vector_3D():
        x(0.), y(0.), z(0.) { }
    vector_3D(float x, float y, float z):
        x(x), y(y), z(z) { }
};

void vector_3D_copy(vector_3D &init, vector_3D &final){
    final.x = init.x;
    final.y = init.y;
    final.z = init.z;
}


struct Atom
{
    vector_3D impulse;
    vector_3D coordinate_0;
    vector_3D coordinate;
    int index;
    unsigned int collisions_laser;
    unsigned int collisions_chaos;
    unsigned int photons_inside;
    unsigned int photons_tmp;
    float radius_2; //квадрат радиуса
    bool tmp;
    Atom(vector_3D impulse,  vector_3D coordinate_0, vector_3D coordinate, 
        unsigned int photon_absorption, unsigned int collisions_laser, unsigned int collisions_chaos, // photon absorption - change cell
        unsigned int photons_inside, unsigned int photons_tmp, float radius_2, bool tmp) :
            impulse{ impulse }, coordinate_0{ coordinate_0 }, coordinate{ coordinate }, index{ index }, collisions_laser{collisions_laser}, collisions_chaos{ collisions_chaos },
            photons_inside{ photons_inside }, photons_tmp{ photons_tmp }, radius_2{radius_2}, tmp{ tmp} { }
    Atom():
        impulse{vector_3D()}, coordinate_0{ vector_3D()}, coordinate{ vector_3D() }, index{ -1 }, collisions_laser{0}, collisions_chaos{ 0 },
        photons_inside{ 0 }, photons_tmp{ 0 }, radius_2{0.0576}, tmp{ false} { } //TODO посчитать квадрат радиуса + проверить численные константы
    void print(){
        std::cout << coordinate.x << ' ' << coordinate.y << ' ' << coordinate.z << ' ' << photons_inside << ' ' << photons_tmp << ' ' << tmp << endl;
    }
    void write(std::ofstream& out){
        out << coordinate_0.x << ' ' << coordinate_0.y << ' ' << coordinate_0.z << ' ';
        out << impulse.x << ' ' << impulse.y << ' ' << impulse.z << ' ';
        out << coordinate.x << ' ' << coordinate.y << ' ' << coordinate.z << ' ';
        out << collisions_laser << ' ' << collisions_chaos << ' ' ;
        out << photons_inside << ' ' << photons_tmp << ' ' << tmp << ' ' << index << std::endl;
    }
    void move(){
        coordinate.x += impulse.x;
        coordinate.y += impulse.y;
        coordinate.z += impulse.z;
    }
};

struct Photon
{
    vector_3D* emission;
    vector_3D* dir;
    Photon(vector_3D* emission, vector_3D* dir): emission{emission}, dir{dir} {}
    Photon(): emission{new vector_3D()}, dir{new vector_3D()} {} 
    void print(){
        std::cout << "photon " << dir->x << ' ' << dir->y << ' ' << dir->z << ' ' << emission->x << ' ' << emission->y << ' ' << emission->z << endl;
    }
    ~Photon(){
        delete emission;
        delete dir;
    }
};



//Шарик, в котором содержится много атомов, чтобы при поиске атома, в который попадет фотон, не перебирать все атомы
// (ищем шары, которые пересекает луч, и только в них перебираем атомы)
struct Ball_cell
{
    vector_3D coordinate_centre;
    Atom* elements;
    unsigned int number_of_elements;
    float c_radius_2;
    Ball_cell():
        coordinate_centre{vector_3D()}, elements{new Atom[cell_size]}, number_of_elements{ 0 }, c_radius_2{ 0 } {}
    void print(){ 
        std::cout << coordinate_centre.x << ' ' << coordinate_centre.y << ' ' << coordinate_centre.z << ' ' << c_radius_2 << ' ' << number_of_elements << endl;
    }
    ~Ball_cell() {
        //print(); 
        //delete[] coordinate_centre; 
        delete[] elements;
    }
};

void write(std::string file, Ball_cell* ball_cells, unsigned const int ball_cells_size, Atom* out_of_cloud, unsigned int* size_out_of_cloud){
    std::ofstream out;
    out.open(file);
    if (out.is_open()) {
        for (unsigned int i = 0; i < ball_cells_size; i++){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }
                ball_cells[i].elements[j].write(out);
            }
        }

        for (unsigned int i = 0; i < *size_out_of_cloud; i++){
            if (out_of_cloud[i].tmp == false){
                std::cout << 127 << endl;
            }
            out_of_cloud[i].write(out);
        }
    }
    out.close();
    return;
}



//Расстояние между точками и точкой и прямой в квадрате
float distance_between_points_2(vector_3D* point_1, vector_3D* point_2){
    //return pow((pow((point_1[0] - point_2[0]), 2) + pow((point_1[1] - point_2[1]), 2) + pow((point_1[2] - point_2[2]), 2)), 0.5);
    //return (point_1[0] - point_2[0])*(point_1[0] - point_2[0]) + (point_1[1] - point_2[1])*(point_1[1] - point_2[1]) + (point_1[2] - point_2[2])*(point_1[2] - point_2[2]);
    return (point_1->x - point_2->x)*(point_1->x - point_2->x) + (point_1->y - point_2->y)*(point_1->y - point_2->y) + (point_1->z - point_2->z)*(point_1->z - point_2->z);
}

float distance_between_line_and_point_2(vector_3D* dir, vector_3D* p_1, vector_3D* p_2){ //TODO проверить формулу
    // return pow(pow((dir[1]*(p_1[2] - p_2[2])-dir[2]*(p_1[1]-p_2[1])), 2)
    //                  + pow((dir[0]*(p_1[2] - p_2[2])-dir[2]*(p_1[0]-p_2[0])), 2)
    //                  + pow((dir[0]*(p_1[1] - p_2[1])-dir[1]*(p_1[0]-p_2[0])), 2), 0.5)/pow(pow(dir[0], 2) + pow(dir[1], 2) + pow(dir[2], 2), 0.5);
    // return ((dir[1]*(p_1[2] - p_2[2])-dir[2]*(p_1[1]-p_2[1])) * (dir[1]*(p_1[2] - p_2[2])-dir[2]*(p_1[1]-p_2[1]))
    //                 + (dir[0]*(p_1[2] - p_2[2])-dir[2]*(p_1[0]-p_2[0])) * (dir[0]*(p_1[2] - p_2[2])-dir[2]*(p_1[0]-p_2[0])) 
    //                 + ((dir[0]*(p_1[1] - p_2[1])-dir[1]*(p_1[0]-p_2[0]))) * ((dir[0]*(p_1[1] - p_2[1])-dir[1]*(p_1[0]-p_2[0])))); // /(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    return ((dir->y*(p_1->z - p_2->z)-dir->z*(p_1->y-p_2->y)) * (dir->y*(p_1->z - p_2->z)-dir->z*(p_1->y-p_2->y))
                    + (dir->x*(p_1->z - p_2->z)-dir->z*(p_1->x-p_2->x)) * (dir->x*(p_1->z - p_2->z)-dir->z*(p_1->x-p_2->x)) 
                    + ((dir->x*(p_1->y - p_2->y)-dir->y*(p_1->x-p_2->x))) * ((dir->x*(p_1->y - p_2->y)-dir->y*(p_1->x-p_2->x)))); // /(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
    
}

float scalar_product(float* a, float *b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

vector_3D generate_random_point(float* random_array, unsigned long long random_array_size, unsigned long long* counter){

    while (true){
        float x = ( random_array[*counter] - 0.5 )* 2;
        (*counter)++;

        if (*counter >= random_array_size){
            *counter = fabs(rand()%(random_array_size-5));
        }


        float y = ( random_array[*counter] - 0.5 )* 2;
        (*counter)++;

        if (*counter >= random_array_size){
            *counter = fabs(rand()%(random_array_size-5));
        }

        float z = ( random_array[*counter] - 0.5 )* 2;
        (*counter)++;

        if (*counter >= random_array_size){
            *counter = fabs(rand()%(random_array_size-5));
        }

        if ( x * x + y * y + z * z > 0.5 && x * x + y * y + z * z < 1){

            return vector_3D(x/sqrt(x * x + y * y + z * z),
                             y/sqrt(x * x + y * y + z * z), z/sqrt(x * x + y * y + z * z));
        }
    }
}


//Столкновение фотона и атома
void collision(Photon* photon, Atom* atom){
    atom->photons_tmp += 1;
    if (photon->dir->x == 1){ 
        atom->collisions_laser += 1;
    }
    else{
        atom->collisions_chaos += 1;
    }

    atom->impulse.x += photon->dir->x * impulse_coefficient;
    atom->impulse.y += photon->dir->y * impulse_coefficient;
    atom->impulse.z += photon->dir->z * impulse_coefficient;

}


//Испускание фотона (закон сохранения импульса)
void emission(Photon* photon, Atom* atom){
    atom->photons_tmp -= 1;
    
    atom->impulse.x -= photon->dir->x * impulse_coefficient;
    atom->impulse.y -= photon->dir->y * impulse_coefficient;
    atom->impulse.z -= photon->dir->z * impulse_coefficient;
    
}


//Если фотон вылетел из облака
void fly_away(Photon* photon, unsigned int number_of_photons, unsigned long long* photon_out, int const screen_size, int const screen_resolution, vector< vector<int> >* screen, int const screen_len){
    (*photon_out) += number_of_photons;
    if (photon->dir->x == 1. ){
        if (photon->emission->y + screen_size < 0 || photon->emission->y + screen_size > 2 * screen_size){
            photon->print();
            cout << 147 << endl;
            
        }
        if (photon->emission->z + screen_size < 0 || photon->emission->z + screen_size > 2 * screen_size){
            photon->print();
            cout << 153 << endl;
            
        }

        (*screen)[static_cast<long long unsigned int>((photon->emission->y + screen_size) * screen_resolution)][static_cast< long long unsigned int>((photon->emission->z + screen_size) * screen_resolution)] += number_of_photons;
    }

}


//Ищем ближайший атом, куда прилетит фотон
Atom* find_nearest_atom(Photon* photon, Atom* atom_departure, Ball_cell* ball_cells, unsigned const int ball_cells_size){

    Atom* nearest_atom = nullptr;
    for (unsigned int i = 0; i < ball_cells_size; i++){ //TODO size_t
        if(ball_cells[i].number_of_elements == 0){
            continue;
        }

        if ((photon->dir->x*(ball_cells[i].coordinate_centre.x - photon->emission->x) + //TODO разность массивов
            photon->dir->y*(ball_cells[i].coordinate_centre.y - photon->emission->y) + 
            photon->dir->z*(ball_cells[i].coordinate_centre.z - photon->emission->z)) < 0
            and distance_between_points_2(&(atom_departure->coordinate), &(ball_cells[i].coordinate_centre)) > 3*ball_cells[i].c_radius_2){
                continue;
        }
        if (nearest_atom != nullptr){
            if (distance_between_points_2(&(nearest_atom->coordinate), &(atom_departure->coordinate)) + 2 * ball_cells[i].c_radius_2 <\
                distance_between_points_2(&(atom_departure->coordinate), &(ball_cells[i].coordinate_centre))){ //TODO функция работает неправильно
                    continue; 
            }
        }

        if (distance_between_line_and_point_2(photon->dir, photon->emission, &(ball_cells[i].coordinate_centre)) < ball_cells[i].c_radius_2){
            
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }

                if (&ball_cells[i].elements[j] == atom_departure){
                    continue;
                }

                if (photon->dir->x * (ball_cells[i].elements[j].coordinate.x - photon->emission->x) + 
                    photon->dir->y * (ball_cells[i].elements[j].coordinate.y - photon->emission->y) + 
                    photon->dir->z * (ball_cells[i].elements[j].coordinate.z - photon->emission->z) < 0){
                        continue;
                }

                if (distance_between_line_and_point_2(photon->dir, photon->emission, &(ball_cells[i].elements[j].coordinate)) < ball_cells[i].elements[j].radius_2){  // 0.24 = (sigma_0/pi)**0.5/sigma_r*r
                    if (not nearest_atom){
                        nearest_atom = &ball_cells[i].elements[j];
                    }
                    else if (distance_between_points_2(&(ball_cells[i].elements[j].coordinate), &(atom_departure->coordinate)) < distance_between_points_2(&(nearest_atom->coordinate), &(atom_departure->coordinate))){
                        nearest_atom = &ball_cells[i].elements[j];
                    }
                }
            }
        }
    }
    return nearest_atom;
}

//
void next_collision(Photon* photon, Atom* atom_departure, unsigned int number_of_photons, Ball_cell* ball_cells, unsigned const int ball_cells_size, unsigned long long* photon_out, int const screen_size, int const screen_resolution, vector< vector<int> >* screen, int const screen_len){
    if (number_of_photons == 1){
        Atom* nearest_atom = find_nearest_atom(photon, atom_departure, ball_cells, ball_cells_size);
       
        if (nearest_atom){
            collision(photon, nearest_atom);
            if (nearest_atom->photons_tmp > 1){

                emission(photon, nearest_atom);
                emission(photon, nearest_atom);
                //Photon* new_photon{new Photon(nearest_atom->coordinate, photon->dir)};
                //photon->emission = nearest_atom->coordinate; //TODO оператор присваивания?
                vector_3D_copy(nearest_atom->coordinate, *(photon->emission));
                

                Atom* next_atom = find_nearest_atom(photon, nearest_atom, ball_cells, ball_cells_size);

                if (not next_atom){
                    fly_away(photon, 2, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 0){
                    collision(photon, next_atom);
                    vector_3D_copy(next_atom->coordinate, *(photon->emission));
                    next_collision(photon, next_atom, 1, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 1){
                    vector_3D_copy(next_atom->coordinate, *(photon->emission));
                    emission(photon, next_atom);
                    next_collision(photon, next_atom, 3, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }
            }
            else {
                return;
            }
        }
        else{
            //std::cout << *photon_out << endl;
            fly_away(photon, 1, photon_out, screen_size, screen_resolution, screen, screen_len);
            //std::cout << *photon_out << endl;
            return;
        }
        return;
    }

    else if (number_of_photons >= 2){
        Atom* nearest_atom = find_nearest_atom(photon, atom_departure, ball_cells, ball_cells_size);
        if (nearest_atom != nullptr){
            vector_3D_copy(nearest_atom->coordinate, *(photon->emission));
            if (nearest_atom->photons_tmp == 0){
                collision(photon, nearest_atom);
                next_collision(photon, nearest_atom, number_of_photons - 1, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                return;
            }
            
            if (nearest_atom->photons_tmp == 1){
                
                emission(photon, nearest_atom);
                Atom* next_atom = find_nearest_atom(photon, nearest_atom, ball_cells, ball_cells_size);

                if (not next_atom){
                    fly_away(photon, number_of_photons+1, photon_out,screen_size, screen_resolution, screen, screen_len);
                    return;
                }
                
                else if (next_atom->photons_tmp == 0){
                    collision(photon, next_atom);
                    vector_3D_copy(next_atom->coordinate, *(photon->emission));
                    next_collision(photon, next_atom, number_of_photons, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 1){
                    emission(photon, next_atom);
                    vector_3D_copy(next_atom->coordinate, *(photon->emission));
                    next_collision(photon, next_atom, number_of_photons+2, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }
            }
        }
        else{
            fly_away(photon, number_of_photons, photon_out, screen_size, screen_resolution, screen, screen_len);
            return;
        }
        return;
    }
}

void emit_photon(Photon* photon, Atom* atom, float* random_array, unsigned long long random_array_size, unsigned long long* counter){

//    float phi = random_array[*counter] * 2 * 3.14;
//    (*counter)++;
//
//    if (*counter >= random_array_size){
//        *counter = fabs(rand()%(random_array_size-5));
//    }
//
//
//    float theta = random_array[*counter]*2*3.14;
//    (*counter)++;
//
//    if (*counter >= random_array_size){
//        *counter = fabs(rand()%(random_array_size-5));
//    }
//
//    vector_3D_copy(atom->coordinate, *(photon->emission));
//
//    photon->dir->x = cos(theta) * cos(phi);
//    photon->dir->y = cos(theta) * sin(phi);
//    photon->dir->z = sin(theta);


    vector_3D direction = generate_random_point(random_array, random_array_size, counter);

    vector_3D_copy(atom->coordinate, *(photon->emission));

    photon->dir->x = direction.x;
    photon->dir->y = direction.y;
    photon->dir->z = direction.z;
    
    atom->impulse.x-= photon->dir->x * impulse_coefficient;
    atom->impulse.y -= photon->dir->y * impulse_coefficient;
    atom->impulse.z -= photon->dir->z * impulse_coefficient;
}


//Если атом вылетел за границы своей клетки, то меняем клетку
void change_cell(Atom* atom, Ball_cell* cell_init, Ball_cell* ball_cells, unsigned const int ball_cells_size, unsigned long long* photon_out, Atom* out_of_cloud, unsigned int* out_of_cloud_size, unsigned long long* photon_in){

    if (distance_between_points_2(&(atom->coordinate), &(cell_init->coordinate_centre)) < cell_init->c_radius_2){
        return;
    }

    atom->tmp = false;
    
    bool break_flag = false;
    for (unsigned int i = 0; i < ball_cells_size; i++){
        if (distance_between_points_2(&(atom->coordinate), &(ball_cells[i].coordinate_centre)) < ball_cells[i].c_radius_2){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements + 2; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    //ball_cells[i].elements[j] = *atom; // TODO Оператор присваивания

                    ball_cells[i].elements[j].tmp = true;
                    vector_3D_copy(atom->coordinate_0, ball_cells[i].elements[j].coordinate_0);
                    vector_3D_copy(atom->coordinate, ball_cells[i].elements[j].coordinate);
                    vector_3D_copy(atom->impulse, ball_cells[i].elements[j].impulse);
                    ball_cells[i].elements[j].photons_inside = atom->photons_inside;
                    ball_cells[i].elements[j].photons_tmp = atom->photons_tmp;
                    ball_cells[i].elements[j].radius_2 = atom->radius_2;
                    ball_cells[i].elements[j].index = atom->index;
                    ball_cells[i].elements[j].collisions_chaos = atom->collisions_chaos;
                    ball_cells[i].elements[j].collisions_laser = atom->collisions_laser;

                    break_flag = true;
                    
                    if (j >= ball_cells[i].number_of_elements){
                        ball_cells[i].number_of_elements += 1;
                    }
                    break;
                }
            }
        }
        if (ball_cells[i].number_of_elements >= cell_size - 1){
            std::cout << 406 << endl;
        }

        if (break_flag){
            break;
        }
    }

    if (break_flag == false){
        
        if (atom->photons_tmp == 1 || atom->photons_inside == 1){
            (*photon_out)++;
            atom->photons_tmp = 0;
            atom->photons_inside = 0;
        }
        out_of_cloud[*out_of_cloud_size].tmp = true;
        vector_3D_copy(atom->coordinate_0, out_of_cloud[*out_of_cloud_size].coordinate_0);
        vector_3D_copy(atom->coordinate, out_of_cloud[*out_of_cloud_size].coordinate);
        vector_3D_copy(atom->impulse, out_of_cloud[*out_of_cloud_size].impulse);
        out_of_cloud[*out_of_cloud_size].photons_inside = atom->photons_inside;
        out_of_cloud[*out_of_cloud_size].photons_tmp = atom->photons_tmp;
        out_of_cloud[*out_of_cloud_size].radius_2 = atom->radius_2;
        out_of_cloud[*out_of_cloud_size].index = atom->index;
        out_of_cloud[*out_of_cloud_size].collisions_chaos = atom->collisions_chaos;
        out_of_cloud[*out_of_cloud_size].collisions_laser = atom->collisions_laser;

        (*out_of_cloud_size)++;  

        if (*out_of_cloud_size >= 10000){
            std::cout << 360 << endl;
        }
    }    
}



void clear_cells(Ball_cell* ball_cells, unsigned int const ball_cells_size){ //TODO написать функцию, которая будет очищать клетки
    for (unsigned int i = 0; i < ball_cells_size; i++){
        if (ball_cells[i].number_of_elements == 0){
            continue;
        }
        while (ball_cells[i].elements[ball_cells[i].number_of_elements-1].tmp != 1){
            ball_cells[i].number_of_elements --;
            if (ball_cells[i].number_of_elements == 0){
                break;
            }
            
        }
        if (ball_cells[i].number_of_elements < 0){
            std::cout << 453 << endl;
        }
    }
}

//Problems:
//destructor
// const ball_cells_size

// intensity_coeff - мощность луча
void experiment(float intensity_coeff, Ball_cell* ball_cells, unsigned int const ball_cells_size, Atom* out_of_cloud, unsigned int* size_out_of_cloud, float const concentration){

    int const screen_size = 56;
    //Детектор фотонов. screen находится далеко за облаком. По сути, это то, что увидит камера (тень от атомов).
    // screen это матрица камеры
    int const screen_resolution = 10; 
    int const screen_len = (screen_size*2*screen_resolution) * (screen_size*2*screen_resolution);
    vector <vector <int> > * screen = new vector <vector <int> >;
    screen->resize(2 * screen_size * screen_resolution);
    for (unsigned int i = 0; i < 2*screen_size*screen_resolution; i++){
        (*screen)[i].resize(2 * screen_size * screen_resolution);
        for (int j = 0; j < 2 * screen_size * screen_resolution; j++){
            (*screen)[i][j] = 0;
        }
    }

    unsigned long long photon_in = 0;
    unsigned long long photon_out = 0;
    int number_of_atoms = 0;
    int number_of_excited_atoms = 0;

    //Массив, заполненный рандомными числами от 0 до 1,
    // нужен для генерирования фотонов от лазера (выбираем рандомную точку в квадрате, ткуда они вылетают)
    // и испускания фотонов от атомов в рандомном направлении (функция emit_photon)
    unsigned long long seed = 50000000;
    float* random_array = new float[seed];
    unsigned long long random_array_size = seed;
    unsigned long long pointer = 0; 
    mt19937 rng(seed); 
    for (unsigned i = 0; i < random_array_size; i++){
        random_array[i] = float(rng())/(pow(2, 32) - 1);
    }

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << intensity_coeff;
    std::string s = stream.str();
    stream << std::fixed << std::setprecision(2) << concentration;
    std::string s1 = stream.str();


    string title = "random_gauss_10.04_" + s1 + "_56_" + s; //title
    mkdir(title.c_str());
    //mkdir(title.c_str(), 0700); //LINUX

    Atom trial_atom;
    Photon photon;
    
    // Количество фотонов, которые вылетают из лазера за итерацию
    int photons_per_iter = int(2760 * intensity_coeff * pow(screen_size - 1, 2)/pow(28, 2));

    time_t start, end;
    time(&start);
    for (int time = 0; time < 10001; time++){
        //Фотоны вылетают из лазера.
        for (int i = 0; i < photons_per_iter; i++){ 
            photon_in++;
            if (pointer + 3 > random_array_size){
                pointer = fabs(rand()%(random_array_size-5));
            }

            trial_atom.coordinate.x = -3*screen_size; 
            trial_atom.coordinate.y = (random_array[pointer++] - 0.5) * 2 * (screen_size-1); 
            trial_atom.coordinate.z = (random_array[pointer++] - 0.5) * 2 * (screen_size-1);


            //Все фотоны от лазера летят вдоль оси X
            photon.dir->x = 1;
            photon.dir->y = 0;
            photon.dir->z = 0;
            vector_3D_copy(trial_atom.coordinate, *(photon.emission));
            next_collision(&photon, &trial_atom, 1, ball_cells, ball_cells_size, &photon_out,  screen_size, screen_resolution, screen, screen_len); //TODO убрать trial_atom
        }

        //Проходимся по атомам, которые расположены в хаотичном порядке. Если в прошлой итерации у атома был фотон, излучаем его.
        //mt19937 rng(1);
        for (unsigned int i = 0; i < ball_cells_size; i++){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp and ball_cells[i].elements[j].photons_inside == 1){
                    ball_cells[i].elements[j].photons_inside = 0;
                    //TODO memcpy
                    vector_3D_copy(ball_cells[i].elements[j].coordinate, *(photon.emission));
                    emit_photon(&photon, &ball_cells[i].elements[j], random_array, random_array_size, &pointer);
                    next_collision(&photon, &ball_cells[i].elements[j], 1, ball_cells, ball_cells_size, &photon_out,  screen_size, screen_resolution, screen, screen_len);
                }
            }
        }

        //Двигаем атомы. Если атом вышел из клетки, помещаем его в новую или в массив атомов, которые вылетели из облака
        for (unsigned int i = 0; i < ball_cells_size; i++){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                
                if (ball_cells[i].elements[j].tmp){
                    ball_cells[i].elements[j].move(); 
                    change_cell(&ball_cells[i].elements[j], &ball_cells[i], ball_cells, ball_cells_size, &photon_out, out_of_cloud, size_out_of_cloud, &photon_in);
                }
                if (ball_cells[i].number_of_elements >= cell_size - 2){
                    std::cout<< 100 << endl;
                }
            }
        }
        for (unsigned int i = 0; i < *size_out_of_cloud; i++){
            out_of_cloud[i].move();
        }

        for (unsigned int i = 0; i < ball_cells_size; i++){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp){
                    ball_cells[i].elements[j].photons_inside = ball_cells[i].elements[j].photons_tmp;
                    ball_cells[i].elements[j].photons_tmp = 0; // TODO segmentation fault
                }
            }
        }

        //Очищаем клетки, которые освободились
        clear_cells(ball_cells, ball_cells_size);

        //Считаем атомы и вылетевшие фотоны для самопроверки
        number_of_atoms = 0;
        number_of_excited_atoms = 0;
        for (unsigned int i = 0; i < ball_cells_size; i++){
            for (unsigned int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }
                number_of_atoms++;
                if (ball_cells[i].elements[j].photons_inside == 1) {
                    number_of_excited_atoms++;
                }
                else if (ball_cells[i].elements[j].photons_inside != 0){
                    std::cout << 545 << endl;
                }
            }
        }
        std::cout << intensity_coeff << ' '<< time << ' ' << photon_in << ' ' << photon_out << ' ' << number_of_excited_atoms << ' ' <<  number_of_excited_atoms + photon_out << ' ' << number_of_atoms << ' ' << *size_out_of_cloud << endl;
    
        //Дальше запись в файл  
        if (time % 50 == 0){ //TODO 250
            write(title + "/movement_data_gauss_"+ to_string(time) + ".txt", ball_cells, ball_cells_size, out_of_cloud, size_out_of_cloud);

            std::ofstream out;
            out.open(title + "/screen_data_gauss_"+ to_string(time) + ".txt");
            if (out.is_open()) {
                for (unsigned int i = 0; i < 2 * screen_size * screen_resolution; i++){
                    for (unsigned int j = 0; j < 2 * screen_size * screen_resolution; j++){
                        out << (*screen)[i][j] << ' ' << i << ' ' << j << std::endl;
                    }
                }
            }
            out.close();

        }


    }


    time(&end);
    double seconds = difftime(end, start);
    printf("The time: %f seconds\n", seconds);

    delete[] random_array;
    delete screen; 
    
    std::cout << 538 << endl;
    delete[] ball_cells; 
    delete[] out_of_cloud;
    std::cout << 540 << endl;

    return;
}


//Здесь считываем координаты атомов и клеток из файлов и запускаем эксперимент.
void read_and_exp(float intensity, float const concentration){
    //Массив клеток, в каждой из которого хранится < cell_size атомов
    unsigned const int ball_cells_size = 20844; 
    Ball_cell* ball_cells = new Ball_cell[ball_cells_size];

    //Массив вылетевших атомов
    Atom* out_of_cloud = new Atom[10000]; //see change_cell 360
    unsigned int size_out_of_cloud = 0;
    
    // Считываем с файла координаты клеток
    ifstream fileptr ("d:\\Programs\\c++\\distributions\\cells_56.txt"); //Менять количество клеток вручную
    //ifstream fileptr ("/mnt/d/Programs/c++/distributions/cells_56.txt"); //LINUX
    if (not fileptr) std::cout << "Файл не открыт!!!\n\n"; //TODO try and catch

    for (unsigned int i = 0; i < ball_cells_size; i++){
        fileptr >> ball_cells[i].coordinate_centre.x >> ball_cells[i].coordinate_centre.y >> ball_cells[i].coordinate_centre.z;
        ball_cells[i].c_radius_2 = (sqrt(3) * 6 / 2) * (sqrt(3) * 6 / 2);
        ball_cells[i].number_of_elements = 0;
    }

    fileptr.close();
    std::cout<< ball_cells[0].coordinate_centre.x<<endl;

    // Добавляем в клетки атомы
    unsigned int number_of_atoms_0 = int(concentration * 73846);

    //Массив хаотично расположенных атомов

    Atom trial_atom;
    int index_counter = 0;
    ifstream file ("d:\\Programs\\c++\\distributions\\atoms_4_56.txt"); //TODO менять количество атомов вручную и размер клетки
    //ifstream file ("/mnt/d/Programs/c++/distributions/atoms_4_56.txt"); //LINUX title
    if (not file) cout << "файл не открыт!!!\n\n";
    for (unsigned int i = 0; i < number_of_atoms_0; i++){
        file >> trial_atom.coordinate_0.x >> trial_atom.coordinate_0.y >> trial_atom.coordinate_0.z;
        // file >> trial_atom.coordinate_0.x >> trial_atom.coordinate_0.y >> trial_atom.coordinate_0.z >> //TODO это разве не функция write?
        //         trial_atom.impulse.x >> trial_atom.impulse.y >> trial_atom.impulse.z >>
        //         trial_atom.coordinate.x >> trial_atom.coordinate.y >> trial_atom.coordinate.z >> 
        //         trial_atom.collisions_laser >> trial_atom.collisions_chaos >>
        //         trial_atom.photons_inside >> trial_atom.photons_tmp >> trial_atom.tmp;
        bool break_flag = true;
        //trial_atom.print();
        for (unsigned int j = 0; j < ball_cells_size; j++){
            if (distance_between_points_2(&trial_atom.coordinate_0, &ball_cells[j].coordinate_centre) <  ball_cells[j].c_radius_2){
                vector_3D_copy(trial_atom.coordinate_0, ball_cells[j].elements[ball_cells[j].number_of_elements].coordinate);
                vector_3D_copy(trial_atom.coordinate_0, ball_cells[j].elements[ball_cells[j].number_of_elements].coordinate_0);
                vector_3D_copy(trial_atom.impulse, ball_cells[j].elements[ball_cells[j].number_of_elements].impulse);
                ball_cells[j].elements[ball_cells[j].number_of_elements].collisions_chaos = trial_atom.collisions_chaos;
                ball_cells[j].elements[ball_cells[j].number_of_elements].collisions_laser = trial_atom.collisions_laser;
                ball_cells[j].elements[ball_cells[j].number_of_elements].index = index_counter++;
                ball_cells[j].elements[ball_cells[j].number_of_elements].tmp = true;
                ball_cells[j].number_of_elements++;
                if (ball_cells[j].number_of_elements > cell_size - 2){
                    cout << 469 << endl;
                }
                break_flag = false;
                break;
            }
        }
        if (break_flag){
            vector_3D_copy(trial_atom.coordinate_0, out_of_cloud[size_out_of_cloud].coordinate);
            vector_3D_copy(trial_atom.coordinate_0, out_of_cloud[size_out_of_cloud].coordinate_0);//TODO pi и sqrt(3)
            vector_3D_copy(trial_atom.impulse, out_of_cloud[size_out_of_cloud].impulse);
            out_of_cloud[size_out_of_cloud].collisions_chaos = trial_atom.collisions_chaos;
            out_of_cloud[size_out_of_cloud].collisions_laser = trial_atom.collisions_laser;
            out_of_cloud[size_out_of_cloud].index = index_counter++;
            out_of_cloud[size_out_of_cloud].tmp = true;
            out_of_cloud[size_out_of_cloud].print();
            size_out_of_cloud++;
        }

    }
    file.close();

    experiment(intensity, ball_cells, ball_cells_size, out_of_cloud, &size_out_of_cloud, concentration);
    return;
}


int main(){
    

    // for (float i = 0.2; i < 1; i += 0.4){
    //    read_and_exp(i);
    // }
    // for (float i = 1; i <= 5; i+=2){
    //     read_and_exp(i);
    // }
    float concentration = 1;
//    for (float i = 1; i <= 1; i += 3){
//        read_and_exp(i, concentration);
//    }
    read_and_exp(4, concentration);
    return 0;
}