
void run_numbers_in_subset(int year, TString RS_data_files)
{

    TFileCollection* fc = new TFileCollection("fc", "fc", RS_data_files);
    TChain* t = new TChain("DecayTree");
    t->AddFileInfoList((TCollection*)fc->GetList());

    ROOT::RDataFrame df(*t);
    std::vector<unsigned int> all_v = df.Take<unsigned int>("runNumber").GetValue();
    sort( all_v.begin(), all_v.end() );

    std::vector<unsigned int> unique_values;

    for(int i = 0; i < size(all_v); i++)
    {
        if(i == 0)
        {
            unique_values.push_back(all_v[i]);
        }
        else
        {
            if(all_v[i] != all_v[i-1])
            {
                unique_values.push_back(all_v[i]);
            }
        }
    }

    std::ofstream file_stream;
    file_stream.open (Form("/panfs/felician/B2Ktautau/workflow/get_run_numbers/201%i/run_numbers_in_subset.txt",year));   

    for(int i = 0; i < size(unique_values); i++)
    {
        file_stream << unique_values[i] << endl;
    }                                           

    return;
}
