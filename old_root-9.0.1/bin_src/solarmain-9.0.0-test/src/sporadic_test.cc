static std::exception_ptr threadExceptionPtr = nullptr;
load_traits(interp, file, headers, Y, 0, headers.size(), ids.size()) 
static void load_trait_set(const int thread_index, const char * phenotype_filename, \
    vector<string> trait_list, double * Y, const int start_index, const int stop_index,\
    const int n_rows){
        const char * errmsg = 0;
        SolarFile * file;
        try{
            file = SolarFile::open("sporadic", phenotype_filename, &errmsg);
            if(errmsg){
                string error_message = "Thread " + to_string(thread_index) + " Error:" + errmsg;
                throw runtime_error(error_message);
            }
            file->start_setup(&errmsg);
            if(errmsg){
                string error_message = "Thread " + to_string(thread_index) + " Error:" + errmsg;
                throw runtime_error(error_message);
            }
            for(int i = start_index; i < stop_index; i++){
                file->setup(trait_list[i],&errmsg);
                if(errmsg){
                    string error_message = "Thread " + to_string(thread_index) + " Error:" + errmsg;
                    throw runtime_error(error_message);
                }                
            }

            const int batch_size = stop_index - start_index;
            int row = 0;
            char ** file_data;
  
            while (0 != (file_data = file->get (&errmsg))){
                if(errmsg){
                    string error_message = "Thread " + to_string(thread_index) + " Error:" + errmsg;
                    throw runtime_error(error_message);
                }                
        
                for(int col = 0; col < batch_size; col++){
                    if(strlen(file_data[col]) != 0){
                        Y[col*n_rows + row] = strtod(file_data[col], NULL);
                    }else{
                        Y[col*n_rows + row] = nan("");
                    }
                }
                row++;
            }                            
        }catch(...)
        {
            threadExceptionPtr = current_exception();
        }
}
static void sporadic_normalize_set(const int thread_index, const int start_index, const int stop_index,  vector<string> ids, vector<string> covariates, unordered_map<string, vector<double> > covariate_map,
                             const char * phenotype_filename, double * Y, double * residuals, const int n_covariates, const int n_rows){
    

    const char * errmsg = 0;
    for(int col_index = start_index; col_index < stop_index; col_index++){
        vector<string> id_list;
        double * current_Y = Y + col_index*n_rows;
        for(int row = 0; row < ids.size(); row++){
            double value = current_Y[row];
            if((value == value && (covariate_map[id].size() != 0)) || (value == value && n_covariates == 0)) {
                id_list.push_back(id);
            }
        }
        Eigen::VectorXd residual;
        if(n_covariates != 0){
            Eigen::MatrixXd covariate_matrix = create_covariate_matrix(id_list, covariates,
                        covariate_map,  n_covariates);
            Eigen::VectorXd col_Y(id_list.size());
            vector<double> covariate_vector;
            int index = 0;
            for(int row = 0 ; row < ids.size(); row++){
                double value = current_Y[row];
                string id = ids[row];
                if(value == value && covariate_map[id].size() != 0)
                    col_Y(index++) = value;

            }
            Eigen::VectorXd beta = covariate_matrix.colPivHouseholderQr().solve(col_Y);

            residual = col_Y - covariate_matrix*beta;
        }else{
            Eigen::VectorXd col_Y(id_list.size());
            int index = 0;
            for(int row = 0 ; row < ids.size(); row++){
                double value = current_Y[row];
                string id = ids[row];
                if(value == value)
                    col_Y(index++) = value;
            }
            residual = col_Y;            
        }
        vector< pair<double, size_t> > data_in;
        unordered_map<size_t, size_t> output_map;
        vector<string>::iterator bad_id_iter;
        int row = 0;
        for(int index = 0 ; index < ids.size(); index++){
            if(find(id_list.begin(), id_list.end(), ids[index]) != id_list.end()){
                data_in.push_back(pair<double, size_t>(residual(row), row));
                output_map[row++] = index;
            }else{
                residuals[index + col_index*n_rows] = nan("");
            }
        }
        inormalize(data_in, output_map, residuals + col_index*n_rows);                
    }
    
}



extern "C" int sporadicNormalizeCmd(ClientData clientData, Tcl_Interp * interp,
                                    int argc, const char * argv[]){

    int max_threads = std::thread::hardware_concurrency();
    const char * header_filename = 0;
    
    const char * out_filename = 0;
    string class_list_str;
   // auto start = std::chrono::high_resolution_clock::now();
    for(int arg = 1 ; arg < argc; arg++){
        if((!StringCmp(argv[arg], "--list", case_ins) || !StringCmp(argv[arg], "-list", case_ins)) && arg + 1 != argc){
            header_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--out", case_ins) || !StringCmp(argv[arg], "-out", case_ins) || !StringCmp(argv[arg], "-o", case_ins) ||\
                  !StringCmp(argv[arg], "--o", case_ins)) && arg + 1 != argc){
            out_filename = argv[++arg];
        }else if((!StringCmp(argv[arg], "--class", case_ins) || !StringCmp(argv[arg], "-class", case_ins) || !StringCmp(argv[arg], "-c", case_ins) ||\
                  !StringCmp(argv[arg], "--c", case_ins)) && arg + 1 != argc){
            class_list_str = string(argv[++arg]);
        }else if(!StringCmp(argv[arg], "help", case_ins)){
            print_help(interp);
            return TCL_OK;
        }else{
            RESULT_LIT("An invalid argument was entered");
            return TCL_ERROR;
        }
    }
    string phenotype_filename =Phenotypes::filenames();
    if(!phenotype_filename.length()){
    RESULT_LIT("No phenotype file is currently loaded");
    return TCL_ERROR;
    }
    if(!out_filename){
        RESULT_LIT("No output file was specified");
        return TCL_ERROR;
    }
    vector<string> headers;
    if(!header_filename){
        if(Trait::Number_Of() == 0){
            RESULT_LIT("No trait was specified or file containing a list of traits");
            return TCL_ERROR;
        }
        
        headers.push_back(string(Trait::Name(0)));
    }else{
        headers = load_headers(header_filename);
        if(!headers.size()){
            RESULT_LIT("No traits could be read in from file specified");
            return TCL_ERROR;
        }
    vector<string> missing_headers = check_headers_in_phenotype( phenotype_filename.c_str(),  headers);
        if(headers.size() == 0){
       RESULT_LIT("Fields listed in header file were not found in loaded phenotype.");
       return TCL_ERROR;
    }
    if(missing_headers.size() != 0){
       cout << "The following fields listed in the header file were not found in the loaded phenotype file: \n";
       for(vector<string>::iterator field_iter = missing_headers.begin(); field_iter != missing_headers.end(); field_iter++){
        cout << " -" << *field_iter << endl;
        }
    }
    }

    
    Covariate * c;
    vector<string> covariate_terms;
    vector<string> ids;
    unordered_map<string, vector<double> > covariate_map;
    int n_covariates = 0;
    for (int i = 0;( c = Covariate::index(i)); i++)
    {
        CovariateTerm * cov_term;
        
        for(cov_term = c->terms(); cov_term; cov_term = cov_term->next){
            bool found = false;
            
            for(vector<string>::iterator cov_iter = covariate_terms.begin(); cov_iter != covariate_terms.end(); cov_iter++){
                if(!StringCmp(cov_term->name, cov_iter->c_str(), case_ins)){
                    found = true;
                    break;
                }
            }
            if(!found){
                covariate_terms.push_back(string(cov_term->name));
            }
        }
        n_covariates++;
    }


    
    
    if(class_list_str.length() != 0){
        vector<int> class_list;
        vector<int> lines_included;
       // char * token = strdup(class_list_str.c_str());
       // token = strtok((char*)class_list_str.c_str(), &token, ",");
      //  int class_number;// = strtol(class_list_str.c_str(), &token, 10);
        //class_list.push_back(class_number);
        int string_index = 0;
        string str_number;
        while(string_index != class_list_str.length()){
            if(class_list_str[string_index] != ','){
                str_number += class_list_str[string_index];
            }else{
                class_list.push_back(atoi(str_number.c_str()));
                str_number.clear();
            }
            string_index++;
        }
        if(str_number.length() != 0){
           class_list.push_back(atoi(str_number.c_str()));
        }
        /*
        while(class_number = strtol(token, &token, 10)){
            //class_number = strtol(token, &token, 10);
            class_list.push_back(class_number);
        }*/
        ofstream output_stream(out_filename);
        output_stream << "ID,Class";
        for(vector<string>::iterator header_iter = headers.begin(); header_iter != headers.end(); header_iter++){
            output_stream << "," << *header_iter;
        }
        output_stream << "\n";
        for(vector<int>::iterator class_iter = class_list.begin(); class_iter != class_list.end(); class_iter++){
            ids.clear();
            covariate_map.clear();
            lines_included.clear();
            if(load_ids_and_covariates_per_class(interp, phenotype_filename.c_str(),  covariate_terms, \
                                       ids,  covariate_map, lines_included, *class_iter) == TCL_ERROR) return TCL_ERROR;
            double * residuals = new double[ids.size()*headers.size()];
            
            double * Y = new double[ids.size()*headers.size()];
            if(headers.size() < max_threads){
                const char * errmsg = 0;
                SolarFile * file = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                if(errmsg){
                    RESULT_LIT(errmsg);
                    return TCL_ERROR;
                }
                if(load_traits_per_inclusion(interp, file, headers, lines_included, Y, 0, headers.size(), ids.size()) == TCL_ERROR){
                    return TCL_ERROR;
                }
                for(int column = 0; column < headers.size(); column++){
                    sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                       phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                }
                
            }else{
#pragma omp parallel
                {
                    int n_threads = omp_get_num_threads();
                        const char * errmsg = 0;
                     int thread_idx = omp_get_thread_num();
                        SolarFile * file  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                        int batch_size = ceil(double(headers.size())/n_threads);
                        int end = (1 +thread_idx)*batch_size;
                        if(thread_idx + 1 == n_threads)
                            end = headers.size();
                        load_traits_per_inclusion(interp, file, headers, lines_included, Y + thread_idx*batch_size*ids.size(), thread_idx*batch_size,  end, ids.size());
                   delete file;
           }
                    
                   
                    
#pragma omp parallel for
                    for(int column = 0; column < headers.size(); column++){
                        sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                           phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                    }
                    
                
 
            }
            delete [] Y;
            for(int row = 0; row < ids.size(); row++){
                output_stream << ids[row] << "," << *class_iter;
                for(int column = 0; column < headers.size(); column++){
                    double value = residuals[column*ids.size() + row];
                    if(value == value)
                        output_stream << "," << value;
                    else
                        output_stream << ",";
                }
                output_stream << endl;
            }
            
            delete [] residuals;
 
        }
        output_stream.close();
        
    }else{
        if(load_ids_and_covariates(interp, phenotype_filename.c_str(),  covariate_terms, \
                                   ids,  covariate_map) == TCL_ERROR) return TCL_ERROR;
        
        /*
         for(int index = 0; index < covariate_terms.size(); index++){
         if(!StringCmp(covariate_terms[index].c_str(), "SEX", case_ins)){
         vector<double> covariate_row;
         for(vector<string>::iterator id_iter = ids.begin(); id_iter != ids.end(); id_iter++){
         string id = *id_iter;
         if(covariate_map[id].size() == 0)
         continue;
         if(covariate_map[id][index] == 2.0){
         covariate_map[id][index]  = 1.0;
         }else{
         covariate_map[id][index] = 0.0;
         }
         }
         break;
         }
         }*/
        double * residuals = new double[ids.size()*headers.size()];
        
        double * Y = new double[ids.size()*headers.size()];
        if(headers.size() < max_threads){
            const char * errmsg = 0;
            SolarFile * file = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
            if(errmsg){
                RESULT_LIT(errmsg);
                return TCL_ERROR;
            }
            if(load_traits(interp, file, headers, Y, 0, headers.size(), ids.size()) == TCL_ERROR){
                return TCL_ERROR;
            }
            for(int column = 0; column < headers.size(); column++){
                sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                   phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
            }
            
        }else{
         //  SolarFile ** files;
            /*
            const char * errmsgg;
            SolarFile * file_test = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsgg);
            auto start = high_resolution_clock::now();
            load_traits(interp, file_test, headers, Y, 0, headers.size(), ids.size());
            delete file_test;
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<microseconds>(stop - start);
            cout << "Time taken by function: " << duration.count() << " microseconds" << endl;
            start = high_resolution_clock::now(); */
           // vector<thread> threads;
           // omp_set_num_threads(7);       
#pragma omp parallel
            {
               int n_threads = omp_get_num_threads();
             
               const char * errmsg = 0;
                SolarFile * file;
//#pragma omp critical
                file  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
              //  for(int f = 0; f < n_threads; f++){
                     int thread_idx  = omp_get_thread_num();
                   
//                    const char * errmsg = 0;
                  //  files[f]  = SolarFile::open("sporadic", phenotype_filename.c_str(), &errmsg);
                    int batch_size = ceil(double(headers.size())/n_threads);
                    int end = (1 +thread_idx)*batch_size;
                    if(thread_idx + 1 == n_threads)
                        end = headers.size();
                    load_traits(interp, file, headers, Y + thread_idx*batch_size*ids.size(), thread_idx*batch_size,  end, ids.size());
              //  }
                 delete file;
             }
            // stop = high_resolution_clock::now();
            // duration = duration_cast<microseconds>(stop - start);  
            // cout << "Time taken by function: " << duration.count() << " microseconds" << endl; 
               // delete [] files;
                
#pragma omp parallel for
                for(int column = 0; column < headers.size(); column++){
                    sporadic_normalize(interp,  ids, covariate_terms, covariate_map,
                                       phenotype_filename.c_str(), Y + column*ids.size(), residuals + column*ids.size(),  n_covariates);
                }
                
           // }
          //  delete [] files;
            
            
            
        }
        
        delete [] Y;
        
        ofstream file_out(out_filename);
        
        
        file_out << "ID";
        for(vector<string>::iterator trait_iter = headers.begin(); trait_iter != headers.end(); trait_iter++){
            file_out << "," << *trait_iter;
        }
        file_out << endl;
        for(int row = 0; row < ids.size(); row++){
            file_out << ids[row];
            for(int column = 0; column < headers.size(); column++){
                double value = residuals[column*ids.size() + row];
                if(value == value)
                    file_out << "," << value;
                else
                    file_out << ",";
            }
            file_out << endl;
        }
        file_out.close();
        delete [] residuals;
    }
    

    
  //  auto elapsed = std::chrono::high_resolution_clock::now() - start;
    
  //  auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(elapsed);
   // cout << seconds.count() << endl;
    return TCL_OK;
    
}