extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

#include "math/matrix.h"
#include "exception.h"

#include "bts/cmd.h"
#include "bts/common.h"
#include "bts/file.h"

#include "bts/fibre/base/set.h"
#include "bts/fibre/track/set.h"

using namespace BTS;

std::vector<std::vector<Fibre::Base::Object> > recreate_objects(
        const std::vector<const char*>& elem_props) {
    
    std::vector<std::vector<Fibre::Base::Object> > objects(10);
    
    for (size_t size_i = 0; size_i < 10; ++size_i) {
        objects[size_i];
        for (size_t value_i = 0; value_i < 10; ++value_i) {
            objects[size_i].push_back(Fibre::Base::Object(size_i, 3, elem_props));
            objects[size_i][value_i].set(value_i + 1001);
        }
    }
    
    return objects;
    
}

std::vector<std::map<const char*, double> > recreate_prop_row_maps(
        const std::vector<const char*>& elem_props) {
    
    std::vector<std::map<const char*, double> > prop_row_maps(10);
    
    for (size_t value_i = 0; value_i < 10; ++value_i) {
        for (size_t key_i = 0; key_i < elem_props.size(); ++key_i)
            prop_row_maps[value_i][elem_props[key_i]] = (double) (value_i + 1) * 0.01
                    + (double) (key_i + 1) * 0.1 + 0.001;
    }
    
    return prop_row_maps;
    
}

std::vector<std::vector<double> > recreate_prop_rows(const std::vector<const char*>& elem_props) {
    
    std::vector<std::vector<double> > prop_rows(10);
    
    for (size_t value_i = 0; value_i < 10; ++value_i) {
        for (size_t key_i = 0; key_i < elem_props.size(); ++key_i)
            prop_rows[value_i].push_back(
                    (double) (value_i + 1) * 0.01 + (double) (key_i + 1) * 0.1 + 0.001);
    }
    
    return prop_rows;
    
}

std::vector<std::map<std::string, std::string> > create_ext_prop_row_maps() {
    
    std::vector<std::map<std::string, std::string> > prop_rows(10);
    
    for (size_t value_i = 0; value_i < 10; ++value_i) {
        for (size_t key_i = 1; key_i < 11; ++key_i)
            prop_rows[value_i]["EXT_ELEM_PROP_" + str(key_i)] = str(
                    (double) (value_i + 1) * 0.01 + (double) key_i * 1.1 + 0.001);
    }
    
    return prop_rows;
    
}

void recreate_elem_prop_objects(const std::vector<const char*>& elem_props,
                                std::vector<std::vector<Fibre::Base::Object> >& objects,
                                std::vector<std::map<const char*, double> >& prop_row_maps,
                                std::vector<std::vector<double> >& prop_rows) {
    
    objects = recreate_objects(elem_props);
    prop_row_maps = recreate_prop_row_maps(elem_props);
    prop_rows = recreate_prop_rows(elem_props);
    
}
SET_VERSION_DEFAULT
;
SET_AUTHOR("Thomas G. Close");
SET_COPYRIGHT(NULL);

DESCRIPTION = {
    "Tests a random sequence of Fibre::Base::Set operations against those performed on a Python version.",
    "",
    NULL
};

ARGUMENTS= {

    Argument()
};

OPTIONS= {

    Option()};

EXECUTE {
    
        size_t SET_SIZE = 5;
        size_t OBJECT_SIZE = 5;
        
        const char* PROP_1 = "PROP_1";
        const char* PROP_2 = "PROP_2";
        const char* PROP_3 = "PROP_3";
        const char* PROP_4 = "PROP_4";
        const char* PROP_5 = "PROP_5";
//        const char* PROP_6 = "PROP_6";
//        const char* PROP_7 = "PROP_7";
//        const char* PROP_8 = "PROP_8";
//        const char* PROP_9 = "PROP_9";
//        const char* PROP_10 = "PROP_10";
        
        const char* ELEM_PROP_1 = "ELEM_PROP_1";
        const char* ELEM_PROP_2 = "ELEM_PROP_2";
        const char* ELEM_PROP_3 = "ELEM_PROP_3";
        const char* ELEM_PROP_4 = "ELEM_PROP_4";
        const char* ELEM_PROP_5 = "ELEM_PROP_5";
//        const char* ELEM_PROP_6 = "ELEM_PROP_6";
//        const char* ELEM_PROP_7 = "ELEM_PROP_7";
//        const char* ELEM_PROP_8 = "ELEM_PROP_8";
//        const char* ELEM_PROP_9 = "ELEM_PROP_9";
//        const char* ELEM_PROP_10 = "ELEM_PROP_10";
        
        std::string EXT_ELEM_PROP_1 = "EXT_ELEM_PROP_1";
        std::string EXT_ELEM_PROP_2 = "EXT_ELEM_PROP_2";
        std::string EXT_ELEM_PROP_3 = "EXT_ELEM_PROP_3";
        std::string EXT_ELEM_PROP_4 = "EXT_ELEM_PROP_4";
        std::string EXT_ELEM_PROP_5 = "EXT_ELEM_PROP_5";
        std::string EXT_ELEM_PROP_6 = "EXT_ELEM_PROP_6";
        std::string EXT_ELEM_PROP_7 = "EXT_ELEM_PROP_7";
        std::string EXT_ELEM_PROP_8 = "EXT_ELEM_PROP_8";
        std::string EXT_ELEM_PROP_9 = "EXT_ELEM_PROP_9";
        std::string EXT_ELEM_PROP_10 = "EXT_ELEM_PROP_10";
        
        std::vector<const char*> props;
        
        props.push_back(PROP_1);
        props.push_back(PROP_2);
        props.push_back(PROP_3);
        props.push_back(PROP_4);
        props.push_back(PROP_5);
        
        std::vector<const char*> elem_props;
        
        elem_props.push_back(ELEM_PROP_1);
        elem_props.push_back(ELEM_PROP_2);
        elem_props.push_back(ELEM_PROP_3);
        elem_props.push_back(ELEM_PROP_4);
        elem_props.push_back(ELEM_PROP_5);
        
        std::vector<std::vector<Fibre::Base::Object> > objects = recreate_objects(elem_props);
        std::vector<std::vector<double> > prop_rows = recreate_prop_rows(elem_props);
        std::vector<std::map<const char*, double> > prop_row_maps = recreate_prop_row_maps(
                elem_props);
        
        std::vector<std::map<std::string, std::string> > ext_prop_row_maps =
                create_ext_prop_row_maps();
        
        //Iniate set
        Fibre::Base::Set<Fibre::Base::Object> set(SET_SIZE, OBJECT_SIZE, 3, props, elem_props);
        
        for (size_t elem_i = 0; elem_i < SET_SIZE; ++elem_i) {
            Fibre::Base::Object ob = set[elem_i];
            for (size_t elem_elem_i = 0; elem_elem_i < OBJECT_SIZE; ++elem_elem_i) {
                Coord t = ob[elem_elem_i];
                for (size_t elem_elem_elem_i = 0; elem_elem_elem_i < 3; ++elem_elem_elem_i)
                    t[elem_elem_elem_i] = elem_i * 100 + elem_elem_i * 10 + elem_elem_elem_i + 1;
            }
            
            for (size_t prop_i = 0; prop_i < 5; ++prop_i)
                ob.prop(prop_i) = (double) elem_i * 100.0 + (double) (prop_i + 1) * .1 + 0.001;
        }
        
        for (size_t prop_i = 0; prop_i < 5; ++prop_i)
            set.prop(prop_i) = (double) (prop_i + 1) * 1.1 + 0.001;
        
        for (size_t eprop_i = 1; eprop_i < 6; ++eprop_i) {
            set.add_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i), "NaN");
            for (size_t row_i = 0; row_i < set.size(); ++row_i)
                set.set_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i),
                        "Ooogle " + str(eprop_i) + "-" + str(row_i), row_i);
        }
        
        //Iniate set
        Fibre::Base::Set<Fibre::Base::Object> set2(SET_SIZE, OBJECT_SIZE, 3, props, elem_props);
        
        for (size_t elem_i = 0; elem_i < SET_SIZE; ++elem_i) {
            Fibre::Base::Object ob = set2[elem_i];
            for (size_t elem_elem_i = 0; elem_elem_i < OBJECT_SIZE; ++elem_elem_i) {
                Coord t = ob[elem_elem_i];
                for (size_t elem_elem_elem_i = 0; elem_elem_elem_i < 3; ++elem_elem_elem_i)
                    t[elem_elem_elem_i] = -(double) (elem_i * 100 + elem_elem_i * 10
                                                     + elem_elem_elem_i + 1);
            }
            
            for (size_t prop_i = 0; prop_i < 5; ++prop_i)
                ob.prop(prop_i) = -(double) elem_i * 100.0 - (double) (prop_i + 1) * .1 - 0.001;
        }
        
        for (size_t prop_i = 0; prop_i < 5; ++prop_i)
            set2.prop(prop_i) = -(double) (prop_i + 1) * 1.1 - 0.001;
        
        for (size_t eprop_i = 6; eprop_i < 11; ++eprop_i) {
            set2.add_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i), "NaN");
            for (size_t row_i = 0; row_i < set2.size(); ++row_i)
                set2.set_extend_elem_prop("EXT_ELEM_PROP_" + str(eprop_i),
                        "Boogle " + str(eprop_i) + "-" + str(row_i), row_i);
        }
        
        std::ofstream out;
        
        out.open("/home/tclose/Code/Tractography/python/fibre/base/test/output/cpp_output.txt");
        
//
//  //====================================================================================================================
//  // Test functions start here
//  //====================================================================================================================
//
//  out << "\n\n-------------------------------------\n            INITIAL SET             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "\n\n-------------------------------------\n            SECOND SET             \n-------------------------------------\n\n" << set2 << "\n\n";
//  set.insert(objects[5][8],1);
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  out << "set.elem_prop_index(ELEM_PROP_4): " << set.elem_prop_index(ELEM_PROP_4) << std::endl;
//  set.add_elem_prop(ELEM_PROP_2,669);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_7): " << set.has_extend_elem_prop(EXT_ELEM_PROP_7) << std::endl;
//  set.add_prop(PROP_6,984);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.free_elem_degree();
//  set.remove_prop(PROP_3);
//  out << "\n\n-------------------------------------\n            5 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.add_elem_prop(ELEM_PROP_4,463);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.copy_extend_elem_prop_row(set2,3,2);
//  set.copy_extend_elem_prop_row(set2,1,1);
//  out << "set.size(): " << set.size() << std::endl;
//  set.add_prop(PROP_9,419);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.prop(PROP_9): " << set.prop(PROP_9) << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.add_elem_prop(ELEM_PROP_5,661);
//  out << "\n\n-------------------------------------\n            10 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.push_back(objects[9][0], ext_prop_row_maps[3]);
//  out << "set.prop_index(PROP_2): " << set.prop_index(PROP_2) << std::endl;
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.erase(1);
//  set.add_elem_prop(ELEM_PROP_4,815);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.prop_index(PROP_2): " << set.prop_index(PROP_2) << std::endl;
//  set.add_prop(PROP_7,434);
//  set.insert(objects[5][8],1);
//  out << "\n\n-------------------------------------\n            15 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  set.free_elem_degree();
//  set.set_extend_elem_prop_row(ext_prop_row_maps[1],4);
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[3],4);
//  set.push_back(objects[5][6], ext_prop_row_maps[7]);
//  set.erase(2);
//  out << "\n\n-------------------------------------\n            20 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.insert(objects[1][5],5);
//  set.push_back(objects[6][1]);
//  out << "set.elem_prop_index(ELEM_PROP_4): " << set.elem_prop_index(ELEM_PROP_4) << std::endl;
//  out << "set.elem_prop_index(ELEM_PROP_4): " << set.elem_prop_index(ELEM_PROP_4) << std::endl;
//  set.resize(1,574,0);
//  set.push_back(objects[0][2]);
//  set.elem_resize(7,511);
//  out << "\n\n-------------------------------------\n            25 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.push_back(objects[7][8]);
//  out << "set.prop(PROP_4): " << set.prop(PROP_4) << std::endl;
//  set.get_extend_elem_prop_row(1);
//  set.elem_resize(8,886);
//  set.add_extend_elem_props(set2);
//  set.add_elem_prop(ELEM_PROP_4,734);
//  out << "\n\n-------------------------------------\n            30 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.insert(objects[8][4],1);
//  set.remove_elem_prop(ELEM_PROP_3);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.insert(objects[8][5],1);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.set_extend_elem_prop(EXT_ELEM_PROP_3,"733",0);
//  out << "\n\n-------------------------------------\n            35 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.erase(2);
//  set.erase(0);
//  set.erase(1);
//  out << "set.prop(PROP_7): " << set.prop(PROP_7) << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_10,"-99");
//  out << "\n\n-------------------------------------\n            40 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_1);
//  set.free_elem_degree();
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_9,0): " << set.get_extend_elem_prop(EXT_ELEM_PROP_9,0) << std::endl;
//  out << "set.has_elem_prop(ELEM_PROP_9): " << set.has_elem_prop(ELEM_PROP_9) << std::endl;
//  set.elem_resize(8,975);
//  set.remove_elem_prop(ELEM_PROP_5);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_6);
//  out << "\n\n-------------------------------------\n            45 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  set.remove_elem_prop(ELEM_PROP_1);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.elem_prop_index(ELEM_PROP_4): " << set.elem_prop_index(ELEM_PROP_4) << std::endl;
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  set.prop_key(6);
//  set.push_back(objects[8][7]);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_2);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_8);
//  out << "\n\n-------------------------------------\n            50 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.add_elem_prop(ELEM_PROP_4,776);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.copy_extend_elem_prop_row(set2,1,1);
//  set.erase(2);
//  set.push_back(objects[8][8], ext_prop_row_maps[2]);
//  set.elem_resize(2,480);
//  out << "\n\n-------------------------------------\n            55 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.add_extend_elem_prop(EXT_ELEM_PROP_8,"-99");
//  out << "set.size(): " << set.size() << std::endl;
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.insert(objects[2][4],0);
//  out << "set.has_elem_prop(ELEM_PROP_7): " << set.has_elem_prop(ELEM_PROP_7) << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[0],0);
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.push_back(objects[2][2]);
//  out << "\n\n-------------------------------------\n            60 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.copy_extend_elem_prop_row(set2,0,2);
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_3,3): " << set.get_extend_elem_prop(EXT_ELEM_PROP_3,3) << std::endl;
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.add_extend_elem_prop(EXT_ELEM_PROP_4,"-99");
//  set.remove_elem_prop(ELEM_PROP_2);
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  set.prop_key(2);
//  set.remove_elem_prop(ELEM_PROP_4);
//  out << "\n\n-------------------------------------\n            65 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[0],1);
//  set.prop_key(2);
//  out << "set.prop_index(PROP_4): " << set.prop_index(PROP_4) << std::endl;
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_3): " << set.has_extend_elem_prop(EXT_ELEM_PROP_3) << std::endl;
//  set.free_elem_degree();
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_7);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.free_elem_degree();
//  out << "\n\n-------------------------------------\n            70 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_8): " << set.has_extend_elem_prop(EXT_ELEM_PROP_8) << std::endl;
//  set.add_extend_elem_props(set2);
//  set.insert(objects[6][7],2);
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  set.copy_extend_elem_prop_row(set2,3,5);
//  set.copy_extend_elem_prop_row(set2,4,2);
//  set.get_extend_elem_prop_row(4);
//  out << "\n\n-------------------------------------\n            75 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.add_prop(PROP_8,140);
//  set.set_extend_elem_prop(EXT_ELEM_PROP_4,"138",2);
//  set.get_extend_elem_prop_row(3);
//  set.get_extend_elem_prop_row(2);
//  set.set_extend_elem_prop_row(ext_prop_row_maps[5],5);
//  out << "\n\n-------------------------------------\n            80 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.set_extend_elem_prop(EXT_ELEM_PROP_5,"667",2);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_2,"-99");
//  set.insert(objects[1][2],2);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_8,"-99");
//  set.remove_prop(PROP_5);
//  out << "\n\n-------------------------------------\n            85 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_6);
//  set.insert(objects[8][9],3);
//  set.resize(1,277,5);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.add_extend_elem_prop(EXT_ELEM_PROP_6,"-99");
//  set.prop_key(0);
//  out << "\n\n-------------------------------------\n            90 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.erase(0);
//  set.prop_key(3);
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.elem_resize(8,926);
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_9);
//  set.remove_prop(PROP_4);
//  out << "\n\n-------------------------------------\n            95 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_5);
//  set.add_prop(PROP_5,902);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.add_prop(PROP_10,182);
//  set.push_back(objects[8][0], ext_prop_row_maps[2]);
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  out << "set.prop_index(PROP_9): " << set.prop_index(PROP_9) << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.add_prop(PROP_4,779);
//  out << "\n\n-------------------------------------\n            100 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_2): " << set.has_extend_elem_prop(EXT_ELEM_PROP_2) << std::endl;
//  set.add_prop(PROP_3,22);
//  set.resize(3,137,0);
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_10): " << set.has_extend_elem_prop(EXT_ELEM_PROP_10) << std::endl;
//  set.set_extend_elem_prop(EXT_ELEM_PROP_8,"308",0);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_10,2): " << set.get_extend_elem_prop(EXT_ELEM_PROP_10,2) << std::endl;
//  set.erase(1);
//  set.erase(1);
//  out << "\n\n-------------------------------------\n            105 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_3,0): " << set.get_extend_elem_prop(EXT_ELEM_PROP_3,0) << std::endl;
//  set.push_back(objects[8][5], ext_prop_row_maps[4]);
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.prop(PROP_1): " << set.prop(PROP_1) << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.get_extend_elem_prop_row(0);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.elem_resize(9,179);
//  out << "\n\n-------------------------------------\n            110 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_10,1): " << set.get_extend_elem_prop(EXT_ELEM_PROP_10,1) << std::endl;
//  set.push_back(objects[9][1], ext_prop_row_maps[3]);
//  set.erase(0);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.prop(PROP_3): " << set.prop(PROP_3) << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.push_back(objects[9][8]);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.remove_prop(PROP_2);
//  set.get_extend_elem_prop_row(0);
//  out << "\n\n-------------------------------------\n            115 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.add_extend_elem_prop(EXT_ELEM_PROP_5,"-99");
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.prop_key(2);
//  set.prop_key(0);
//  set.erase(0);
//  set.copy_extend_elem_prop_row(set2,4,1);
//  out << "\n\n-------------------------------------\n            120 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.free_elem_degree();
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.remove_prop(PROP_8);
//  out << "set.prop_index(PROP_10): " << set.prop_index(PROP_10) << std::endl;
//  set.free_elem_degree();
//  out << "set.size(): " << set.size() << std::endl;
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_7);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_1,"-99");
//  out << "\n\n-------------------------------------\n            125 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.elem_resize(4,873);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_3);
//  out << "set.size(): " << set.size() << std::endl;
//  set.add_extend_elem_props(set2);
//  set.free_elem_degree();
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_2,0): " << set.get_extend_elem_prop(EXT_ELEM_PROP_2,0) << std::endl;
//  set.insert(objects[5][3],1);
//  out << "\n\n-------------------------------------\n            130 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.push_back(objects[4][7], ext_prop_row_maps[3]);
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.push_back(objects[4][0], ext_prop_row_maps[7]);
//  set.push_back(objects[6][6], ext_prop_row_maps[2]);
//  set.elem_resize(7,192);
//  set.insert(objects[7][9],1);
//  out << "\n\n-------------------------------------\n            135 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.prop(PROP_3): " << set.prop(PROP_3) << std::endl;
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  set.insert(objects[7][1],0);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_6,"-99");
//  out << "set.has_elem_prop(ELEM_PROP_3): " << set.has_elem_prop(ELEM_PROP_3) << std::endl;
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.free_elem_degree();
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.prop(PROP_5): " << set.prop(PROP_5) << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_4): " << set.has_extend_elem_prop(EXT_ELEM_PROP_4) << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.elem_resize(6,44);
//  out << "\n\n-------------------------------------\n            140 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.get_extend_elem_prop_row(2);
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[0],6);
//  set.add_extend_elem_props(set2);
//  set.resize(9,561,0);
//  out << "set.prop_index(PROP_3): " << set.prop_index(PROP_3) << std::endl;
//  set.push_back(objects[6][8]);
//  out << "\n\n-------------------------------------\n            145 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.insert(objects[6][0],5);
//  set.add_extend_elem_props(set2);
//  set.push_back(objects[6][7]);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.prop_key(7);
//  set.push_back(objects[6][1], ext_prop_row_maps[4]);
//  out << "\n\n-------------------------------------\n            150 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  out << "set.has_elem_prop(ELEM_PROP_10): " << set.has_elem_prop(ELEM_PROP_10) << std::endl;
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  set.resize(5,969,0);
//  out << "set.has_prop(PROP_8): " << set.has_prop(PROP_8) << std::endl;
//  set.elem_resize(2,367);
//  set.add_prop(PROP_2,692);
//  set.add_prop(PROP_8,147);
//  set.resize(8,245,0);
//  out << "\n\n-------------------------------------\n            155 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_5,3): " << set.get_extend_elem_prop(EXT_ELEM_PROP_5,3) << std::endl;
//  set.prop_key(6);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_8);
//  set.add_extend_elem_props(set2);
//  set.get_extend_elem_prop_row(5);
//  set.set_extend_elem_prop(EXT_ELEM_PROP_2,"640",0);
//  out << "\n\n-------------------------------------\n            160 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.free_elem_degree();
//  set.get_extend_elem_prop_row(0);
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_9): " << set.has_extend_elem_prop(EXT_ELEM_PROP_9) << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_2);
//  out << "set.prop(PROP_9): " << set.prop(PROP_9) << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  set.free_elem_degree();
//  set.freeze_elem_degree();
//  out << "\n\n-------------------------------------\n            165 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_6);
//  set.remove_prop(PROP_9);
//  set.push_back(objects[2][6]);
//  set.remove_prop(PROP_8);
//  out << "\n\n-------------------------------------\n            170 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.resize(1,823,0);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  set.elem_resize(1,218);
//  set.prop_key(4);
//  out << "set.size(): " << set.size() << std::endl;
//  set.copy_extend_elem_prop_row(set2,4,0);
//  out << "set.prop_index(PROP_3): " << set.prop_index(PROP_3) << std::endl;
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.erase(0);
//  out << "\n\n-------------------------------------\n            175 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_7);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.elem_resize(2,672);
//  set.prop_key(0);
//  set.prop_key(2);
//  set.elem_resize(4,946);
//  out << "\n\n-------------------------------------\n            180 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.prop_key(7);
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_4);
//  out << "set.size(): " << set.size() << std::endl;
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  set.add_prop(PROP_8,859);
//  set.add_prop(PROP_9,835);
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  set.add_extend_elem_prop(EXT_ELEM_PROP_10,"-99");
//  out << "\n\n-------------------------------------\n            185 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.has_prop(PROP_9): " << set.has_prop(PROP_9) << std::endl;
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  out << "set.has_prop(PROP_2): " << set.has_prop(PROP_2) << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_8): " << set.has_extend_elem_prop(EXT_ELEM_PROP_8) << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_3): " << set.has_extend_elem_prop(EXT_ELEM_PROP_3) << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.resize(9,576,0);
//  set.insert(objects[4][4],1);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_7,"-99");
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_9,6): " << set.get_extend_elem_prop(EXT_ELEM_PROP_9,6) << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[6],8);
//  out << "\n\n-------------------------------------\n            190 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.remove_prop(PROP_6);
//  set.set_extend_elem_prop(EXT_ELEM_PROP_7,"493",2);
//  set.elem_resize(9,337);
//  set.add_prop(PROP_6,153);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  set.free_elem_degree();
//  out << "\n\n-------------------------------------\n            195 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.prop_index(PROP_3): " << set.prop_index(PROP_3) << std::endl;
//  set.add_extend_elem_props(set2);
//  set.erase(4);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_1,"-99");
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_5,4): " << set.get_extend_elem_prop(EXT_ELEM_PROP_5,4) << std::endl;
//  out << "set.prop_index(PROP_10): " << set.prop_index(PROP_10) << std::endl;
//  set.erase(8);
//  set.set_extend_elem_prop_row(ext_prop_row_maps[0],2);
//  out << "\n\n-------------------------------------\n            200 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  out << "set.size(): " << set.size() << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_2): " << set.has_extend_elem_prop(EXT_ELEM_PROP_2) << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_1): " << set.has_extend_elem_prop(EXT_ELEM_PROP_1) << std::endl;
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_9,3): " << set.get_extend_elem_prop(EXT_ELEM_PROP_9,3) << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  out << "set.has_prop(PROP_2): " << set.has_prop(PROP_2) << std::endl;
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  out << "set.has_prop(PROP_1): " << set.has_prop(PROP_1) << std::endl;
//  set.copy_extend_elem_prop_row(set2,0,1);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  set.get_extend_elem_prop_row(6);
//  out << "set.num_props(): " << set.num_props() << std::endl;
//  out << "set.elem_prop_keys(): " << set.elem_prop_keys() << std::endl;
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.add_extend_elem_props(set2);
//  out << "\n\n-------------------------------------\n            205 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.elem_resize(2,89);
//  set.push_back(objects[2][6], ext_prop_row_maps[3]);
//  set.elem_resize(5,359);
//  set.free_elem_degree();
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.remove_extend_elem_prop(EXT_ELEM_PROP_10);
//  out << "\n\n-------------------------------------\n            210 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  out << "set.has_extend_elem_prop(EXT_ELEM_PROP_1): " << set.has_extend_elem_prop(EXT_ELEM_PROP_1) << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  out << "set.prop_index(PROP_7): " << set.prop_index(PROP_7) << std::endl;
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.insert(objects[6][7],5);
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.set_extend_elem_prop_row(ext_prop_row_maps[3],0);
//  set.add_extend_elem_prop(EXT_ELEM_PROP_6,"-99");
//  set.remove_prop(PROP_4);
//  out << "set.num_extend_elem_props(): " << set.num_extend_elem_props() << std::endl;
//  set.add_extend_elem_props(set2);
//  out << "\n\n-------------------------------------\n            215 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.num_elem_props(): " << set.num_elem_props() << std::endl;
//  set.add_extend_elem_prop(EXT_ELEM_PROP_5,"-99");
//  set.remove_prop(PROP_7);
//  set.add_extend_elem_props(set2);
//  set.elem_resize(5,955);
//  set.resize(4,369,0);
//  out << "\n\n-------------------------------------\n            220 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.add_prop(PROP_4,364);
//  set.freeze_elem_degree();
//  elem_props = set.get_elem_prop_header(); recreate_elem_prop_objects(elem_props,objects,prop_row_maps,prop_rows);
//  set.add_prop(PROP_7,495);
//  set.free_elem_degree();
//  out << "set.size(): " << set.size() << std::endl;
//  set.push_back(objects[6][6]);
//  out << "\n\n-------------------------------------\n            225 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  out << "set.size(): " << set.size() << std::endl;
//  set.free_elem_degree();
//  out << "set.vsize(): " << set.vsize() << std::endl;
//  out << "set.get_extend_elem_prop(EXT_ELEM_PROP_5,3): " << set.get_extend_elem_prop(EXT_ELEM_PROP_5,3) << std::endl;
//  set.append_extend_elem_props(set2);
//  out << "set.prop_keys(): " << set.prop_keys() << std::endl;
//  set.push_back(objects[0][4]);
//  set.add_extend_elem_props(set2);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.get_extend_elem_prop_header(): " << set.get_extend_elem_prop_header() << std::endl;
//  out << "set.prop(PROP_10): " << set.prop(PROP_10) << std::endl;
//  set.copy_extend_elem_prop_row(set2,3,4);
//  out << "\n\n-------------------------------------\n            230 MODFICATIONS             \n-------------------------------------\n\n" << set << "\n\n";
//  set.get_extend_elem_prop_row(1);
//  out << "set.prop_index(PROP_10): " << set.prop_index(PROP_10) << std::endl;
//  out << "set.has_elem_prop(ELEM_PROP_8): " << set.has_elem_prop(ELEM_PROP_8) << std::endl;
//  set.set_extend_elem_prop(EXT_ELEM_PROP_10,"898",1);
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  out << "set.bsize(): " << set.bsize() << std::endl;
//  set.get_extend_elem_prop_row(1);
        
        out << std::endl << std::endl << "-------------------------------------" << std::endl;
        out << "             FINAL PRINT             " << std::endl;
        out << "-------------------------------------\n" << std::endl;
        
        out << set << std::endl;
        
        std::cout << "Ran cpp test successfully" << std::endl;
        
    }
