//
//  display_phi2.cc
//  
//
//  Created by Brian Donohue on 12/22/21.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "solar.h"
using namespace std;
static const char * write_phi2_graph(const char * output_filename, Tcl_Interp * interp){
    
    const char * errmsg = 0;
   
    SolarFile * pedindex = SolarFile::open("display_phi2", "pedindex.out", &errmsg);
    if(errmsg) return errmsg;
    pedindex->start_setup(&errmsg);
    if(errmsg) return errmsg;
    pedindex->setup("ID", &errmsg);
    if(errmsg) return errmsg;
    int n_ibdids = 1;
    char ** file_data;
    
    while (0 != (file_data = pedindex->get (&errmsg))){
        n_ibdids++;
    }

    Matrix * phi2 = 0;
    phi2 = Matrix::find("phi2");
    if(!phi2){
        Solar_Eval(interp, "matrix load phi2.gz phi2");
        phi2 = Matrix::find("phi2");
        if(!phi2){
            return "phi2 matrix could not be loaded";
        }
    }
    ofstream out_stream(output_filename);

    for(int col = 1; col < n_ibdids; col++){
        
        for(int row = 1 ; row < n_ibdids; row++){
            double phi2_value = 0;
            phi2_value = phi2->get(col, row);
            out_stream << phi2_value*phi2_value;
            if(row + 1 != n_ibdids){
                out_stream << ",";
            }else{
                out_stream << "\n";
            }
        }
    }
   
    out_stream.close();
   
    delete pedindex;
    return 0;
    
}
static void print_phi2_help(Tcl_Interp * interp){
    RESULT_LIT("help print_phi2");
}
extern "C" int print_phi2(ClientData clientData, Tcl_Interp * interp,
                            int argc, const char * argv[]){
    const char * output_filename = 0;
    for(int arg = 1 ;arg < argc ; arg++){
        if(!StringCmp(argv[arg], "help", case_ins) || !StringCmp(argv[arg], "-help", case_ins) || !StringCmp(argv[arg], "--help", case_ins)
           || !StringCmp(argv[arg], "h", case_ins) || !StringCmp(argv[arg], "-h", case_ins) || !StringCmp(argv[arg], "--help", case_ins)){
            print_phi2_help(interp);
            return TCL_OK;
        }else if ((!StringCmp(argv[arg], "-o", case_ins) || !StringCmp(argv[arg], "--o", case_ins) || !StringCmp(argv[arg], "-out", case_ins) ||
                  !StringCmp(argv[arg], "--out", case_ins)) && arg + 1 < argc){
            output_filename = argv[++arg];
        }else{
            RESULT_LIT("Invalid argument enter see help");
            return TCL_ERROR;
        }
    }
    if(output_filename == 0){
        RESULT_LIT("No filename has been entered with -o option");
        return TCL_ERROR;
    }
    const bool is_loaded = loadedPed();
    if(!is_loaded){
        RESULT_LIT("No pedigree has been loaded");
        return TCL_ERROR;
    }

    const char * errmsg = 0 ;
   
    errmsg = write_phi2_graph(output_filename, interp);
    if(errmsg){
        RESULT_BUF(errmsg);
        return TCL_ERROR;
    }
    return TCL_OK;
    
}
